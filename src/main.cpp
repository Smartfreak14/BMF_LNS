#include "BMF.hpp"
#include "BMFLocalSearch.hpp"
#include "CSVMatrixLoader.hpp"
#include <iostream>
#include <iomanip>
#include <random>
#include <map>
#include <set>
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <array>
#include <csignal>
#include <atomic>
#include <sstream>
#include <climits>
#include <tuple>

namespace {

std::atomic<bool> g_interrupted(false);

void signal_handler(int signum) {
    if (signum == SIGINT) {
        std::cout << "\n[interrupt] Ctrl+C received, stopping." << std::endl;
        g_interrupted.store(true);
    }
}

inline bool is_interrupted() { return g_interrupted.load(); }
inline void reset_interrupt() { g_interrupted.store(false); }

}

// phase 1 = greedy ; phase 2 = boucle WLS + LNS partial.
struct CSVTestResult {
    std::string filename;
    std::string method;
    int k = 0;
    int init_errors = 0;
    int phase1_errors = 0;
    double phase1_last_improve_ms = 0.0;
    int phase2_errors = 0;
    double phase2_last_improve_ms = 0.0;
    int final_errors = 0;
    int iterations = 0;
    double total_time = 0;
    long long wls_total_flips = 0;
    double    wls_total_time_ms = 0.0;
    int verif_errors = -1;
    bool success = false;
};

void save_matrix_csv(const Matrix& M, const std::string& filepath) {
    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "save_matrix_csv: cannot open " << filepath << std::endl;
        return;
    }
    for (int i = 0; i < M.rows; i++) {
        for (int j = 0; j < M.cols; j++) {
            file << M(i, j);
            if (j < M.cols - 1) file << ",";
        }
        file << "\n";
    }
}

void print_matrix_stats(const Matrix& M, const std::string& name) {
    int total = M.rows * M.cols;
    int ones = 0, zeros = 0, missing = 0;
    for (int i = 0; i < M.rows; i++) {
        for (int j = 0; j < M.cols; j++) {
            if (M(i, j) == 1) ones++;
            else if (M(i, j) == 0) zeros++;
            else if (M(i, j) == -1) missing++;
        }
    }
    std::cout << "[stats] " << name << " " << M.rows << "x" << M.cols
              << " | ones=" << ones << " (" << std::fixed << std::setprecision(2)
              << (100.0 * ones / total) << "%)"
              << " zeros=" << zeros
              << " missing=" << missing << std::endl;
}

int verify_with_python(const Matrix& M, const Matrix& A, const Matrix& B,
                       const std::string& method_name, const std::string& csv_filename) {
    std::string temp_dir = "/tmp/bmf_verify";
    std::filesystem::create_directories(temp_dir);

    std::string base_name = csv_filename.empty() ? "matrix"
                                                 : std::filesystem::path(csv_filename).stem().string();
    std::string suffix = "_" + method_name;
    std::replace(suffix.begin(), suffix.end(), ' ', '_');
    std::replace(suffix.begin(), suffix.end(), '+', '_');

    std::string A_path = temp_dir + "/" + base_name + suffix + "_A.csv";
    std::string B_path = temp_dir + "/" + base_name + suffix + "_B.csv";
    std::string M_path = temp_dir + "/" + base_name + "_M.csv";

    save_matrix_csv(A, A_path);
    save_matrix_csv(B, B_path);
    save_matrix_csv(M, M_path);

    std::string cmd = "python3 ../src/verif.py \"" + A_path + "\" \"" + B_path + "\" \"" + M_path + "\" 2>&1";

    std::array<char, 256> buf;
    std::string output;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return -1;
    while (fgets(buf.data(), buf.size(), pipe) != nullptr) output += buf.data();
    pclose(pipe);

    const std::string prefix = "Reconstruction error =";
    auto pos = output.find(prefix);
    if (pos == std::string::npos) return -1;
    try {
        return std::stoi(output.substr(pos + prefix.length()));
    } catch (...) {
        return -1;
    }
}

std::tuple<Matrix, int, std::string> load_csv_matrix(const std::string& csv_file, int k_value) {
    std::string filename = csv_file;
    if (filename.find('/') == std::string::npos) {
        filename = "../data/" + filename;
        if (filename.find(".csv") == std::string::npos) filename += ".csv";
    }

    std::cout << "[load] " << filename << std::endl;
    Matrix M = CSVMatrixLoader::loadFromCSV(filename);
    if (M.rows == 0 || M.cols == 0) {
        std::cerr << "load_csv_matrix: cannot load matrix" << std::endl;
        return {Matrix(0, 0), 0, filename};
    }
    print_matrix_stats(M, "M");

    int k = k_value;
    if (k <= 0) {
        k = std::max(5, std::min(50, static_cast<int>(std::log2(std::min(M.rows, M.cols)) * 2)));
        std::cout << "[load] auto k=" << k << std::endl;
    } else {
        std::cout << "[load] k=" << k << std::endl;
    }
    return {M, k, filename};
}

std::pair<std::vector<int>, std::vector<int>> select_neighborhood_topk_random(
    const Matrix& M, const Matrix& A, const Matrix& B, int count, std::mt19937& rng) {
    int m = M.rows, n = M.cols;
    int nr = std::min(count, m), nc = std::min(count, n);

    Matrix computed = A.multiply(B);
    std::vector<std::pair<int, int>> re(m), ce(n);
    for (int i = 0; i < m; i++) {
        int e = 0;
        for (int j = 0; j < n; j++) if (M(i, j) != -1 && computed(i, j) != M(i, j)) e++;
        re[i] = {e, i};
    }
    for (int j = 0; j < n; j++) {
        int e = 0;
        for (int i = 0; i < m; i++) if (M(i, j) != -1 && computed(i, j) != M(i, j)) e++;
        ce[j] = {e, j};
    }
    std::sort(re.begin(), re.end(), std::greater<>());
    std::sort(ce.begin(), ce.end(), std::greater<>());

    int top_r = std::max(1, static_cast<int>(nr * 0.7));
    int rnd_r = nr - top_r;
    int top_c = std::max(1, static_cast<int>(nc * 0.7));
    int rnd_c = nc - top_c;

    std::set<int> sr, sc;
    for (int i = 0; i < top_r && i < m; i++) sr.insert(re[i].second);
    for (int j = 0; j < top_c && j < n; j++) sc.insert(ce[j].second);

    std::vector<int> rem_r, rem_c;
    for (int i = top_r; i < m; i++) rem_r.push_back(re[i].second);
    for (int j = top_c; j < n; j++) rem_c.push_back(ce[j].second);
    std::shuffle(rem_r.begin(), rem_r.end(), rng);
    std::shuffle(rem_c.begin(), rem_c.end(), rng);
    for (int i = 0; i < rnd_r && i < (int)rem_r.size(); i++) sr.insert(rem_r[i]);
    for (int j = 0; j < rnd_c && j < (int)rem_c.size(); j++) sc.insert(rem_c[j]);

    return {std::vector<int>(sr.begin(), sr.end()), std::vector<int>(sc.begin(), sc.end())};
}


namespace {

void setup_signals() {
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
}

auto make_remaining_fn(std::chrono::high_resolution_clock::time_point start, double timeout_ms) {
    return [start, timeout_ms]() -> double {
        if (timeout_ms <= 0) return std::numeric_limits<double>::max();
        auto now = std::chrono::high_resolution_clock::now();
        return timeout_ms - std::chrono::duration<double, std::milli>(now - start).count();
    };
}

inline double elapsed_ms(std::chrono::high_resolution_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t0).count();
}

void announce(const std::string& header) {
    std::cout << "\n=== " << header << " ===" << std::endl;
}

void summarize(const CSVTestResult& r) {
    std::cout << "[" << r.method << "] errors " << r.init_errors
              << " -> p1=" << r.phase1_errors
              << " -> p2=" << r.phase2_errors
              << " | last_improve p1=" << std::fixed << std::setprecision(1) << r.phase1_last_improve_ms << "ms"
              << " p2=" << r.phase2_last_improve_ms << "ms"
              << " | iters=" << r.iterations
              << " | total=" << r.total_time << " ms" << std::endl;
}

} 

// opt4-wls = GI + (WLS + LNS partial)*
CSVTestResult test_opt4_wls(const std::string& csv_file, int k_value, double timeout_seconds) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.method = "opt4-wls";
    result.k = k_value;

    setup_signals();
    announce("opt4-wls : GREEDY + (WLS + LNS partial)*");

    auto [M, k, filename] = load_csv_matrix(csv_file, k_value);
    if (M.rows == 0) return result;

    auto t0 = std::chrono::high_resolution_clock::now();
    auto remaining = make_remaining_fn(t0, timeout_seconds * 1000.0);

    BMFLocalSearch ls(k, M, 42);
    ls.initialize_greedy();
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    double t_after_greedy = elapsed_ms(t0);
    std::cout << "[greedy] init_errors=" << init_errors << std::endl;

    int phase1_errors = init_errors;
    double phase1_last_improve = t_after_greedy;
    int best = init_errors;
    Matrix A_best = ls.A, B_best = ls.B;
    Matrix A_cur = ls.A, B_cur = ls.B;

    int loops = 0, lns_iters = 0;
    std::mt19937 rng(42);
    int nh = 10;
    int nh_max = std::min(100, std::min(M.rows, M.cols));

    BMF bmf(k, M);
    bool phase1_done = false;
    int phase2_errors = phase1_errors;
    double phase2_last_improve = phase1_last_improve;

    while (best > 0 && !is_interrupted() && remaining() > 1000) {
        loops++;

        BMFLocalSearch ls_w(k, M, 42 + loops);
        ls_w.A = A_cur; ls_w.B = B_cur;
        double t_wls_start = elapsed_ms(t0);
        double wls_budget = std::max(1000.0, remaining() * 0.5);
        LocalSearchResult wr = ls_w.solve_weighted(
            200000, true, &g_interrupted, wls_budget,
            -1.0, 1.0, 1.0);

        result.wls_total_flips   += wr.iterations;
        result.wls_total_time_ms += wr.total_time;

        if (wr.final_errors < best) {
            best = wr.final_errors; A_best = wr.A_solution; B_best = wr.B_solution;
            double last_improve = t_wls_start + wr.last_improvement_ms;
            if (!phase1_done) { phase1_errors = best; phase1_last_improve = last_improve; }
            else              { phase2_errors = best; phase2_last_improve = last_improve; }
        } else if (!phase1_done) {
            phase1_errors = best;
            phase1_last_improve = t_wls_start + wr.last_improvement_ms;
        }
        if (!phase1_done) {
            phase2_errors = phase1_errors;
            phase2_last_improve = phase1_last_improve;
            phase1_done = true;
        }

        A_cur = wr.A_solution; B_cur = wr.B_solution;
        if (best == 0 || is_interrupted() || remaining() < 1000) break;

        auto [rows, cols] = select_neighborhood_topk_random(M, A_cur, B_cur, nh, rng);
        int c = bmf.lns_step_partial(A_cur, B_cur, rows, cols);
        lns_iters++;

        if (c >= 0 && c < best) {
            best = c; A_best = A_cur; B_best = B_cur;
            phase2_errors = best;
            phase2_last_improve = elapsed_ms(t0);
        } else {
            nh = std::min(nh_max, nh + 10);
        }
        A_cur = A_best; B_cur = B_best;
    }

    double t_total = elapsed_ms(t0);
    result.k = k;
    result.init_errors = init_errors;
    result.phase1_errors = phase1_errors;
    result.phase1_last_improve_ms = phase1_last_improve;
    result.phase2_errors = phase2_errors;
    result.phase2_last_improve_ms = phase2_last_improve;
    result.final_errors = best;
    result.iterations = lns_iters;
    result.total_time = t_total;
    result.success = true;
    summarize(result);

    if (result.final_errors == 0) {
        result.verif_errors = verify_with_python(M, A_best, B_best, "OPT4_WLS", filename);
    }
    return result;
}

namespace {

static std::string verif_label(const CSVTestResult& r) {
    if (r.verif_errors < 0) return "n/a";
    if (r.verif_errors == 0) return "0";
    return "MISMATCH(" + std::to_string(r.verif_errors) + ")";
}

void print_result_table(const CSVTestResult& r) {
    const int W = 152;
    auto fps = [](const CSVTestResult& x) -> double {
        double s = x.wls_total_time_ms / 1000.0;
        return (s > 0.0) ? (double)x.wls_total_flips / s : 0.0;
    };
    std::cout << "\n" << std::string(W, '=') << "\n"
              << "RESULT\n"
              << std::string(W, '=') << "\n"
              << std::left
              << std::setw(13) << "Method"
              << std::setw(22) << "Dataset"
              << std::right
              << std::setw(7)  << "k"
              << std::setw(8)  << "init"
              << std::setw(11) << "phase1 err"
              << std::setw(12) << "phase1 t_ms"
              << std::setw(11) << "phase2 err"
              << std::setw(12) << "phase2 t_ms"
              << std::setw(9)  << "iters"
              << std::setw(12) << "total t_ms"
              << std::setw(12) << "wls flips"
              << std::setw(13) << "wls flips/s"
              << std::setw(13) << "verif" << "\n"
              << std::string(W, '-') << "\n";
    if (r.success) {
        std::cout << std::left
                  << std::setw(13) << r.method
                  << std::setw(22) << r.filename
                  << std::right
                  << std::setw(7)  << r.k
                  << std::setw(8)  << r.init_errors
                  << std::setw(11) << r.phase1_errors
                  << std::setw(12) << std::fixed << std::setprecision(1) << r.phase1_last_improve_ms
                  << std::setw(11) << r.phase2_errors
                  << std::setw(12) << std::fixed << std::setprecision(1) << r.phase2_last_improve_ms
                  << std::setw(9)  << r.iterations
                  << std::setw(12) << std::fixed << std::setprecision(1) << r.total_time
                  << std::setw(12) << r.wls_total_flips
                  << std::setw(13) << std::fixed << std::setprecision(0) << fps(r)
                  << std::setw(13) << verif_label(r) << "\n";
    }
    std::cout << std::string(W, '=') << "\n";
}

void write_csv_result(const std::string& path, const CSVTestResult& r) {
    std::ofstream f(path);
    if (!f.is_open()) {
        std::cerr << "write_csv_result: cannot open " << path << std::endl;
        return;
    }
    f << "method,dataset,k,init_errors,"
         "phase1_errors,phase1_last_improve_ms,"
         "phase2_errors,phase2_last_improve_ms,"
         "final_errors,iterations,time_ms,"
         "wls_total_flips,wls_total_time_ms,wls_flips_per_sec,"
         "verif_errors\n";
    f << std::fixed << std::setprecision(1);
    if (r.success) {
        double s = r.wls_total_time_ms / 1000.0;
        double fps = (s > 0.0) ? (double)r.wls_total_flips / s : 0.0;
        f << r.method << "," << r.filename << "," << r.k << ","
          << r.init_errors << ","
          << r.phase1_errors << "," << r.phase1_last_improve_ms << ","
          << r.phase2_errors << "," << r.phase2_last_improve_ms << ","
          << r.final_errors << "," << r.iterations << ","
          << r.total_time << ","
          << r.wls_total_flips << "," << r.wls_total_time_ms << "," << fps << ","
          << r.verif_errors << "\n";
    }
    std::cout << "[csv] wrote " << path << std::endl;
}

void print_usage() {
    std::cout
        << "Usage:\n"
        << "  ./rbac_main <file.csv> <k> [timeout_seconds] [--csv-out <path>]\n"
        << "  ./rbac_main --help\n\n"
        << "opt4-wls = GI + (WLS + LNS partial)*. timeout=0 : pas de timeout, Ctrl+C stoppe.\n";
}

} // anonymous namespace

int main(int argc, char* argv[]) {
    if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
        print_usage();
        return (argc < 2) ? 1 : 0;
    }

    std::string csv_out;
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--csv-out" && i + 1 < argc) { csv_out = argv[++i]; continue; }
        positional.push_back(a);
    }
    if (positional.size() < 2) {
        std::cerr << "Erreur : il faut au moins <file.csv> <k>." << std::endl;
        print_usage();
        return 1;
    }

    std::string file = positional[0];
    int k = std::stoi(positional[1]);
    double timeout = (positional.size() >= 3) ? std::stod(positional[2]) : 0;

    CSVTestResult r = test_opt4_wls(file, k, timeout);
    if (!r.success) return 1;

    print_result_table(r);
    if (!csv_out.empty()) write_csv_result(csv_out, r);
    return 0;
}
