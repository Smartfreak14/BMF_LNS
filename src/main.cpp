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
#include <unistd.h>

namespace {

std::atomic<bool> g_interrupted(false);

void signal_handler(int signum) {
    // SIGTERM = kill par SLURM sur timeout ; SIGINT = Ctrl+C.
    if (signum == SIGINT || signum == SIGTERM) {
        std::cout << "\n[interrupt] signal " << signum
                  << " received, stopping (CSV partiel sera ecrit)." << std::endl;
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
    int m = 0;                  // rows de M
    int n = 0;                  // cols de M
    unsigned seed = 42;         // graine RNG
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
    double time_to_zero_ms = -1.0;  // -1 si jamais atteint 0
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

std::tuple<int, int, int> verify_full(const Matrix& M, const Matrix& A, const Matrix& B,
                                       const std::string& method_name, const std::string& csv_filename) {
    // Un dossier par processus : deux runs du meme dataset sur un noeud ne s'ecrasent pas.
    std::string temp_dir = "/tmp/bmf_verify/" + std::to_string(getpid());
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
    if (!pipe) return {-1, -1, -1};
    while (fgets(buf.data(), buf.size(), pipe) != nullptr) output += buf.data();
    pclose(pipe);
    std::error_code ec;
    std::filesystem::remove_all(temp_dir, ec);

    auto parse_int_after = [&](const std::string& prefix) -> int {
        auto pos = output.find(prefix);
        if (pos == std::string::npos) return -1;
        try { return std::stoi(output.substr(pos + prefix.length())); }
        catch (...) { return -1; }
    };
    int total = parse_int_after("Reconstruction error =");
    int cov0  = parse_int_after("Covered 0-entries (hard violations) =");
    int unc1  = parse_int_after("Uncovered 1-entries (soft violations, err_M) =");
    return {total, cov0, unc1};
}

int verify_with_python(const Matrix& M, const Matrix& A, const Matrix& B,
                       const std::string& method_name, const std::string& csv_filename) {
    auto [total, cov0, unc1] = verify_full(M, A, B, method_name, csv_filename);
    if (cov0 > 0) {
        std::cerr << "  [VERIF] WARNING: from-below VIOLATED : "
                  << cov0 << " covered 0-entries + "
                  << unc1 << " uncovered 1-entries" << std::endl;
    }
    return total;
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
    std::signal(SIGINT,  signal_handler);
    std::signal(SIGTERM, signal_handler);   
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

enum class AblationMode { FULL, NO_LNS, NO_DOUBLE_FLIP, NO_WLS_NO_DOUBLE,
                          NO_WEIGHTING, GREEDY_ONLY };

static const char* ablation_name(AblationMode m) {
    switch (m) {
        case AblationMode::FULL:              return "boolmaxls";
        case AblationMode::NO_LNS:            return "boolmaxls-no-lns";
        case AblationMode::NO_DOUBLE_FLIP:    return "boolmaxls-no-double";
        case AblationMode::NO_WLS_NO_DOUBLE:  return "boolmaxls-no-double-no-lns";
        case AblationMode::NO_WEIGHTING:      return "boolmaxls-no-weighting";
        case AblationMode::GREEDY_ONLY:       return "boolmaxls-greedy";
    }
    return "boolmaxls";
}

void write_csv_result(const std::string& path, const CSVTestResult& r,
                      const std::string& regime = "", const std::string& trace_file = "");

// Facteurs 0/1 en CSV : A = m lignes de l valeurs, B = l lignes de n valeurs, l <= k ;
// les k - l facteurs restants sont vides.
bool load_factors(const std::string& a_path, const std::string& b_path, Matrix& A, Matrix& B) {
    auto read = [](const std::string& path, std::vector<std::vector<int>>& rows) {
        std::ifstream f(path);
        if (!f.is_open()) return false;
        std::string line;
        while (std::getline(f, line)) {
            std::vector<int> row;
            std::stringstream ss(line);
            std::string cell;
            while (std::getline(ss, cell, ',')) row.push_back(cell == "1" ? 1 : 0);
            rows.push_back(row);
        }
        return true;
    };
    std::vector<std::vector<int>> a, b;
    if (!read(a_path, a) || !read(b_path, b)) {
        std::cerr << "Erreur : --init-factors illisible : " << a_path << " / " << b_path << std::endl;
        return false;
    }
    int l = (int)b.size();
    bool ok = (int)a.size() == A.rows && l <= A.cols;
    for (const auto& r : a) ok = ok && (int)r.size() == l;
    for (const auto& r : b) ok = ok && (int)r.size() == B.cols;
    if (!ok) {
        std::cerr << "Erreur : --init-factors : A doit etre " << A.rows << "x l et B l x " << B.cols
                  << " avec l <= " << A.cols << std::endl;
        return false;
    }
    for (int i = 0; i < A.rows; i++)
        for (int f = 0; f < l; f++) A(i, f) = a[i][f];
    for (int f = 0; f < l; f++)
        for (int j = 0; j < B.cols; j++) B(f, j) = b[f][j];
    std::cout << "[init-factors] " << l << " facteurs sur " << A.cols << " depuis " << a_path << std::endl;
    return true;
}

CSVTestResult test_boolmaxls(const std::string& csv_file, int k_value, double timeout_seconds,
                            AblationMode ablation = AblationMode::FULL, unsigned seed = 42,
                            const std::string& trace_path = "",
                            const std::string& csv_out = "", const std::string& regime = "",
                            const std::string& init_a = "", const std::string& init_b = "",
                            double lns_limit_s = -1) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.method = ablation_name(ablation);
    result.k = k_value;
    result.seed = seed;

    
    std::ofstream trace;
    if (!trace_path.empty()) {
        trace.open(trace_path);
        if (trace.is_open()) trace << "time_ms err\n";
    }
    auto trace_write = [&](double t_ms, int err) {
        if (trace.is_open()) trace << std::fixed << std::setprecision(3) << t_ms << " " << err << "\n";
    };

    setup_signals();
    announce(std::string(ablation_name(ablation)) + " : ablation=" +
             ablation_name(ablation) + " seed=" + std::to_string(seed));

    auto [M, k, filename] = load_csv_matrix(csv_file, k_value);
    if (M.rows == 0) return result;
    result.m = M.rows;
    result.n = M.cols;

    auto t0 = std::chrono::high_resolution_clock::now();
    auto remaining = make_remaining_fn(t0, timeout_seconds * 1000.0);

    BMFLocalSearch ls(k, M, seed);
    if (init_a.empty()) ls.initialize_greedy();
    else if (!load_factors(init_a, init_b, ls.A, ls.B)) return result;
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    double t_after_greedy = elapsed_ms(t0);
    std::cout << (init_a.empty() ? "[greedy]" : "[init-factors]") << " init_errors=" << init_errors << std::endl;
    trace_write(t_after_greedy, init_errors);

    int phase1_errors = init_errors;
    double phase1_last_improve = t_after_greedy;
    int best = init_errors;
    Matrix A_best = ls.A, B_best = ls.B;
    Matrix A_cur = ls.A, B_cur = ls.B;
    if (best == 0) result.time_to_zero_ms = t_after_greedy;

    auto checkpoint = [&](const std::string& reason) {
        if (csv_out.empty()) return;
        result.init_errors = init_errors;
        result.phase1_errors = phase1_errors;
        result.phase1_last_improve_ms = phase1_last_improve;
        result.phase2_errors = best;
        result.phase2_last_improve_ms = elapsed_ms(t0);
        result.final_errors = best;
        result.total_time = elapsed_ms(t0);
        result.success = true;
        write_csv_result(csv_out, result, regime, trace_path);
        std::cout << "  [checkpoint] " << reason << " best=" << best
                  << " t=" << result.total_time << "ms -> " << csv_out << std::endl;
    };
    checkpoint("init");

    if (ablation == AblationMode::GREEDY_ONLY) {
        result.init_errors = init_errors;
        result.phase1_errors = phase1_errors;
        result.phase1_last_improve_ms = phase1_last_improve;
        result.phase2_errors = phase1_errors;
        result.phase2_last_improve_ms = phase1_last_improve;
        result.final_errors = best;
        result.iterations = 0;
        result.total_time = elapsed_ms(t0);
        result.success = true;
        result.verif_errors = verify_with_python(M, A_best, B_best, result.method, filename);
        summarize(result);
        return result;
    }

    int loops = 0, lns_iters = 0, lns_aborted = 0;
    double lns_total_ms = 0.0, lns_max_ms = 0.0;
    std::mt19937 rng(seed);
    int nh = 10;
    int nh_max = std::min(100, std::min(M.rows, M.cols));

    BMF bmf(k, M);
    bool phase1_done = false;
    int phase2_errors = phase1_errors;
    double phase2_last_improve = phase1_last_improve;

    // NO_WEIGHTING : ni DF, ni LNS, ni ponderation
    bool enable_double = (ablation != AblationMode::NO_DOUBLE_FLIP &&
                          ablation != AblationMode::NO_WLS_NO_DOUBLE &&
                          ablation != AblationMode::NO_WEIGHTING);
    bool enable_lns = (ablation != AblationMode::NO_LNS &&
                       ablation != AblationMode::NO_WLS_NO_DOUBLE &&
                       ablation != AblationMode::NO_WEIGHTING);
    bool enable_weighting = (ablation != AblationMode::NO_WEIGHTING);

    while (best > 0 && !is_interrupted() && remaining() > 1000) {
        loops++;
        int best_before = best;

        BMFLocalSearch ls_w(k, M, seed + loops);
        ls_w.A = A_cur; ls_w.B = B_cur;
        double t_wls_start = elapsed_ms(t0);
        double wls_budget = std::max(1000.0, remaining() * 0.5);
        LocalSearchResult wr = ls_w.solve_weighted(
            200000, true, &g_interrupted, wls_budget,
            -1.0, 1.0, 1.0, enable_double, enable_weighting);

        result.wls_total_flips   += wr.iterations;
        result.wls_total_time_ms += wr.total_time;

        if (wr.final_errors < best) {
            best = wr.final_errors; A_best = wr.A_solution; B_best = wr.B_solution;
            double last_improve = t_wls_start + wr.last_improvement_ms;
            trace_write(last_improve, best);
            if (best == 0 && result.time_to_zero_ms < 0) result.time_to_zero_ms = last_improve;
            if (!phase1_done) { phase1_errors = best; phase1_last_improve = last_improve; }
            else              { phase2_errors = best; phase2_last_improve = last_improve; }
            checkpoint("wls");
        } else if (!phase1_done) {
            phase1_errors = best;
            phase1_last_improve = t_wls_start + wr.last_improvement_ms;
        }
        if (!phase1_done) {
            phase2_errors = phase1_errors;
            phase2_last_improve = phase1_last_improve;
            phase1_done = true;
        }

        A_cur = std::move(wr.A_solution); B_cur = std::move(wr.B_solution);
        if (best == 0 || is_interrupted() || remaining() < 1000) break;

        // Hill-climb pur : relancer depuis le meme minimum local ne change rien.
        if (!enable_weighting && best >= best_before) break;

        if (!enable_lns) {
            continue;
        }

        auto [rows, cols] = select_neighborhood_topk_random(M, A_cur, B_cur, nh, rng);
        double t_lns = elapsed_ms(t0);
        // --lns-time-limit : appel borne par la limite et par le budget restant.
        int c = (lns_limit_s < 0)
                    ? bmf.lns_step_partial(A_cur, B_cur, rows, cols)
                    : bmf.lns_step_partial_limited(A_cur, B_cur, rows, cols,
                                                   std::min(lns_limit_s * 1000.0, remaining()), &g_interrupted);
        double d_lns = elapsed_ms(t0) - t_lns;
        lns_iters++;
        lns_total_ms += d_lns;
        lns_max_ms = std::max(lns_max_ms, d_lns);
        if (c == -2) {
            lns_aborted++;
            std::cout << "  [PARTIAL] abandonne apres " << static_cast<long long>(d_lns) << " ms" << std::endl;
        }

        if (c >= 0 && c < best) {
            best = c; A_best = A_cur; B_best = B_cur;
            phase2_errors = best;
            phase2_last_improve = elapsed_ms(t0);
            trace_write(phase2_last_improve, best);
            if (best == 0 && result.time_to_zero_ms < 0) result.time_to_zero_ms = phase2_last_improve;
            checkpoint("lns");
        } else {
            nh = std::min(nh_max, nh + 10);
        }
        A_cur = A_best; B_cur = B_best;
    }
    if (lns_iters > 0) {
        std::ostringstream s;
        s << std::fixed << std::setprecision(1) << "[LNS] " << lns_iters << " appels, moyenne "
          << lns_total_ms / lns_iters << " ms, max " << lns_max_ms << " ms, abandonnes " << lns_aborted;
        std::cout << s.str() << std::endl;
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

    result.verif_errors = verify_with_python(M, A_best, B_best, result.method, filename);
    return result;
}


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

void write_csv_result(const std::string& path, const CSVTestResult& r,
                      const std::string& regime, const std::string& trace_file) {
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
         "verif_errors,"
         "m,n,regime,tool,seed,time_to_zero_ms,trace_file\n";
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
          << r.verif_errors << ","
          << r.m << "," << r.n << "," << regime << "," << r.method << ","
          << r.seed << "," << (r.time_to_zero_ms < 0 ? -1.0 : r.time_to_zero_ms) << ","
          << trace_file << "\n";
    }
    std::cout << "[csv] wrote " << path << std::endl;
}

void print_usage() {
    std::cout
        << "Usage:\n"
        << "  ./BoolMaxLS <file.csv> <k> [timeout_seconds] [OPTIONS]\n"
        << "  ./BoolMaxLS --dump-wcnf <file.csv> <k> <output.wcnf> [--init-out <file>] [--max-clauses N]\n"
        << "  ./BoolMaxLS --help\n\n"
        << "BoolMaxLS = GI + (WLS + LNS partial)*.  timeout=0 : pas de timeout.\n\n"
        << "OPTIONS :\n"
        << "  --csv-out <path>         Ecrit une ligne CSV de resultat.\n"
        << "  --trace-out <path>       Ecrit la trace (time_ms err) a chaque amelioration.\n"
        << "  --regime <str>           Etiquette regime dans le CSV (raw/simpluni/simplexi/synth).\n"
        << "  --seed <n>               Graine RNG (defaut 42, multi-seed reproducibilite).\n"
        << "  --init-factors <A> <B>   Depart depuis ces facteurs CSV (au plus k, ex. GreEss) au lieu du glouton.\n"
        << "  --ablation <mode>        Modes ablation (defaut full) :\n"
        << "     full                    GI + WLS + LNS  (BoolMaxLS complet)\n"
        << "     greedy                  GI seul\n"
        << "     no-double-flip          GI + WLS(sans double flip) + LNS\n"
        << "     no-lns                  GI + WLS + DF (sans LNS)\n"
        << "     wls-no-double-no-lns    GI + WLS pondere (sans DF, sans LNS)\n"
        << "     no-weighting            GI + WLS non pondere (hill-climb pur, sans DF, sans LNS)\n"
        << "  --lns-time-limit <s|inf>  Chaque sous-probleme exact du LNS est resolu dans un processus fils,\n"
        << "                           abandonne apres s secondes et au plus tard a la fin du budget (inf :\n"
        << "                           seulement a la fin du budget). Sans cette option, aucune limite.\n"
        << "  --dump-wcnf              Ecrit la formule Partial MaxSAT (WCNF v1).\n"
        << "     --init-out <file>       + affectation gloutonne (A, B, auxiliaires) pour nuwls -init.\n"
        << "     --max-clauses N         au-dela de N clauses : en-tete seul (formule non ecrite).\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
        print_usage();
        return (argc < 2) ? 1 : 0;
    }

    // --dump-wcnf <file.csv> <k> <output.wcnf> [--init-out <file>] [--max-clauses N]
    if (std::string(argv[1]) == "--dump-wcnf") {
        if (argc < 5) {
            std::cerr << "Usage: ./BoolMaxLS --dump-wcnf <file.csv> <k> <output.wcnf>"
                         " [--init-out <file>] [--max-clauses N]\n";
            return 1;
        }
        std::string init_out;
        long long max_clauses = 0;
        for (int i = 5; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--init-out" && i + 1 < argc)         { init_out = argv[++i]; continue; }
            if (a == "--max-clauses" && i + 1 < argc)      { max_clauses = std::stoll(argv[++i]); continue; }
            std::cerr << "Erreur : option inconnue " << a << std::endl;
            return 1;
        }
        auto [M, k_out, filename] = load_csv_matrix(argv[2], std::stoi(argv[3]));
        if (M.rows == 0) return 1;
        BMF bmf(k_out, M);
        if (!bmf.dump_wcnf(argv[4], max_clauses)) return 0;
        if (!init_out.empty()) {
            BMFLocalSearch ls(k_out, M, 42);
            ls.initialize_greedy();
            ls.compute_all_counts();
            std::cout << "[init] greedy errors=" << ls.count_errors() << std::endl;
            bmf.dump_assignment(ls.A, ls.B, init_out);
        }
        return 0;
    }

    std::string csv_out, trace_out, regime, init_a, init_b;
    AblationMode ablation = AblationMode::FULL;
    unsigned seed = 42;
    double lns_limit_s = -1;
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--init-factors" && i + 2 < argc) { init_a = argv[++i]; init_b = argv[++i]; continue; }
        if (a == "--csv-out"   && i + 1 < argc) { csv_out   = argv[++i]; continue; }
        if (a == "--trace-out" && i + 1 < argc) { trace_out = argv[++i]; continue; }
        if (a == "--regime"    && i + 1 < argc) { regime    = argv[++i]; continue; }
        if (a == "--seed" && i + 1 < argc) { seed = (unsigned)std::stoul(argv[++i]); continue; }
        if (a == "--ablation" && i + 1 < argc) {
            std::string v = argv[++i];
            if      (v == "full" || v == "boolmaxls") ablation = AblationMode::FULL;
            else if (v == "no-lns")                            ablation = AblationMode::NO_LNS;
            else if (v == "no-double-flip" || v == "no-double") ablation = AblationMode::NO_DOUBLE_FLIP;
            else if (v == "wls-no-double-no-lns")              ablation = AblationMode::NO_WLS_NO_DOUBLE;
            else if (v == "no-weighting")                       ablation = AblationMode::NO_WEIGHTING;
            else if (v == "greedy" || v == "greedy-only")      ablation = AblationMode::GREEDY_ONLY;
            else { std::cerr << "Erreur : --ablation inconnu : " << v << std::endl; return 1; }
            continue;
        }
        if (a == "--lns-time-limit" && i + 1 < argc) {
            std::string v = argv[++i];
            try { lns_limit_s = (v == "inf") ? std::numeric_limits<double>::infinity() : std::stod(v); }
            catch (...) { lns_limit_s = 0; }
            if (!(lns_limit_s > 0)) { std::cerr << "Erreur : --lns-time-limit attend s > 0 ou inf" << std::endl; return 1; }
            continue;
        }
        // Une option mal ecrite serait sinon prise pour un argument positionnel et ignoree.
        if (a.rfind("--", 0) == 0) { std::cerr << "Erreur : option inconnue " << a << std::endl; return 1; }
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

    CSVTestResult r = test_boolmaxls(file, k, timeout, ablation, seed, trace_out, csv_out, regime,
                                     init_a, init_b, lns_limit_s);
    if (!r.success) return 1;

    print_result_table(r);
    if (!csv_out.empty()) write_csv_result(csv_out, r, regime, trace_out);
    return 0;
}
