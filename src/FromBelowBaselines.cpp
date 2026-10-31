// From-below BMF baselines, stopped after at most k factors:
//   grecond    : GreConD, Algorithm 2 of Belohlavek & Vychodil, JCSS 76 (2010)
//   greess     : GreEss, Algorithms 1-2 of Belohlavek & Trnecka, JCSS 81 (2015)
//   mdlgrecond : MDLGreConD, Algorithm 1 of Makhalova & Trnecka, ADAC 15 (2021)
//
// Unknown cells (--missing):
//   free (default) : 1 in the formal context but not in the set U of entries to cover, so a
//                    factor may cover them without having to, as in BoolMaxLS and verif.py
//   zero / one     : replaced by 0 / 1 before the algorithm runs
// Errors are counted on known cells only, as in verif.py.

#include "CSVMatrixLoader.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

using Word = std::uint64_t;

struct Bits {
    std::vector<Word> w;
    int size = 0;

    Bits() = default;
    explicit Bits(int n, bool full = false) : w((n + 63) / 64, full ? ~Word(0) : 0), size(n) {
        if (full) trim();
    }
    void trim() {
        if (size % 64 && !w.empty()) w.back() &= (Word(1) << (size % 64)) - 1;
    }
    void fill() {
        std::fill(w.begin(), w.end(), ~Word(0));
        trim();
    }
    bool get(int i) const { return (w[i >> 6] >> (i & 63)) & 1; }
    void set(int i) { w[i >> 6] |= Word(1) << (i & 63); }
    int count() const {
        int c = 0;
        for (Word x : w) c += __builtin_popcountll(x);
        return c;
    }
    bool none() const {
        for (Word x : w) if (x) return false;
        return true;
    }
    bool subset_of(const Bits& o) const {
        for (size_t t = 0; t < w.size(); t++) if (w[t] & ~o.w[t]) return false;
        return true;
    }
    bool intersects(const Bits& o) const {
        for (size_t t = 0; t < w.size(); t++) if (w[t] & o.w[t]) return true;
        return false;
    }
    Bits& operator&=(const Bits& o) {
        for (size_t t = 0; t < w.size(); t++) w[t] &= o.w[t];
        return *this;
    }
    Bits& operator|=(const Bits& o) {
        for (size_t t = 0; t < w.size(); t++) w[t] |= o.w[t];
        return *this;
    }
    bool operator==(const Bits& o) const { return w == o.w; }
    void assign_and(const Bits& a, const Bits& b) {
        for (size_t t = 0; t < w.size(); t++) w[t] = a.w[t] & b.w[t];
    }
    template <class F> void for_each(F f) const {
        for (size_t t = 0; t < w.size(); t++)
            for (Word x = w[t]; x; x &= x - 1) f(int(t * 64 + __builtin_ctzll(x)));
    }
};

inline int and_count(const Bits& a, const Bits& b) {
    int c = 0;
    for (size_t t = 0; t < a.w.size(); t++) c += __builtin_popcountll(a.w[t] & b.w[t]);
    return c;
}

// Formal context: objects x attributes, stored by rows and by columns.
struct Context {
    int m = 0, n = 0;
    std::vector<Bits> row, col;

    Context(int m_, int n_) : m(m_), n(n_), row(m_, Bits(n_)), col(n_, Bits(m_)) {}
    void set(int i, int j) {
        row[i].set(j);
        col[j].set(i);
    }
    // C^up : attributes shared by all objects of C
    void up(const Bits& c, Bits& d) const {
        d.fill();
        c.for_each([&](int i) { d &= row[i]; });
    }
    // D^down : objects having all attributes of D
    void down(const Bits& d, Bits& c) const {
        c.fill();
        d.for_each([&](int j) { c &= col[j]; });
    }
    Bits up(const Bits& c) const {
        Bits d(n);
        up(c, d);
        return d;
    }
    Bits down(const Bits& d) const {
        Bits c(m);
        down(d, c);
        return c;
    }
};

// Entries still to be covered, by rows.
struct Uncovered {
    std::vector<Bits> row;
    std::vector<int> cnt;
    long long total = 0;

    explicit Uncovered(const std::vector<Bits>& rows) : row(rows), cnt(rows.size()) {
        for (size_t i = 0; i < rows.size(); i++) total += cnt[i] = rows[i].count();
    }
    long long area(const Bits& c, const Bits& d) const {
        long long s = 0;
        c.for_each([&](int i) { if (cnt[i]) s += and_count(row[i], d); });
        return s;
    }
    long long row_sum(const Bits& c) const {
        long long s = 0;
        c.for_each([&](int i) { s += cnt[i]; });
        return s;
    }
    void remove(const Bits& c, const Bits& d) {
        c.for_each([&](int i) {
            if (!cnt[i]) return;
            int before = cnt[i];
            for (size_t t = 0; t < d.w.size(); t++) row[i].w[t] &= ~d.w[t];
            cnt[i] = row[i].count();
            total -= before - cnt[i];
        });
    }
};

struct Clock {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double limit = 0;
    double elapsed() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }
    bool over() const { return limit > 0 && elapsed() > limit; }
};

struct Factor {
    Bits ext, itt;
};

struct Run {
    std::vector<Factor> factors;
    std::vector<long long> uncovered_known;  // [l] = known 1s uncovered by the first l factors
    std::vector<double> time_s;              // [l] = time when factor l was added ([0] = setup)
    std::string stop = "k_reached";
};

// Common bookkeeping: U drives the algorithm, K tracks known 1s for the error.
struct State {
    Uncovered U, K;
    Run run;
    const Clock& clk;

    State(const std::vector<Bits>& ones, const std::vector<Bits>& known_ones, const Clock& c)
        : U(ones), K(known_ones), clk(c) {}
    void start() {
        run.uncovered_known.push_back(K.total);
        run.time_s.push_back(clk.elapsed());
    }
    void add(const Bits& c, const Bits& d) {
        run.factors.push_back({c, d});
        U.remove(c, d);
        K.remove(c, d);
        run.uncovered_known.push_back(K.total);
        run.time_s.push_back(clk.elapsed());
    }
};

// GreConD: Algorithm 2 of Belohlavek & Vychodil (2010).
void grecond(const Context& I, State& st, int k) {
    st.start();
    Bits a(I.m), b(I.n), best_c(I.m), best_d(I.n);
    while (st.U.total > 0) {
        if ((int)st.run.factors.size() >= k) return;
        if (st.clk.over()) { st.run.stop = "timeout"; return; }
        Bits c(I.m, true), d(I.n);  // D = {} and C = D^down
        long long v = 0;
        for (;;) {
            long long best = v;
            int best_j = -1;
            for (int j = 0; j < I.n; j++) {
                if (d.get(j)) continue;
                a.assign_and(c, I.col[j]);  // (D u {j})^down
                long long ca = a.count();
                if (ca == 0 || ca * I.n <= best || st.U.row_sum(a) <= best) continue;
                I.up(a, b);
                long long val = st.U.area(a, b);
                if (val > best) {
                    best = val;
                    best_j = j;
                    best_c = a;
                    best_d = b;
                }
            }
            if (best_j < 0) break;
            c = best_c;
            d = best_d;
            v = best;
        }
        st.add(c, d);
    }
    st.run.stop = "exact";
}

// E(I) by Lemma 2 of Belohlavek & Trnecka (2015).
Context essential_part(const Context& I) {
    std::vector<Bits> rmask(I.m, Bits(I.n)), cmask(I.n, Bits(I.m));
    std::vector<int> rc(I.m), cc(I.n);
    for (int i = 0; i < I.m; i++) rc[i] = I.row[i].count();
    for (int j = 0; j < I.n; j++) cc[j] = I.col[j].count();
    for (int i = 0; i < I.m; i++)
        for (int i2 = 0; i2 < I.m; i2++)
            if (rc[i2] < rc[i] && I.row[i2].subset_of(I.row[i])) rmask[i] |= I.row[i2];
    for (int j = 0; j < I.n; j++)
        for (int j2 = 0; j2 < I.n; j2++)
            if (cc[j2] < cc[j] && I.col[j2].subset_of(I.col[j])) cmask[j] |= I.col[j2];
    Context E(I.m, I.n);
    for (int i = 0; i < I.m; i++)
        I.row[i].for_each([&](int j) {
            if (!rmask[i].get(j) && !cmask[j].get(i)) E.set(i, j);
        });
    return E;
}

// Interval I_{C,D} of B(I), kept as the context J = I n (D^down x C^up) of Lemma 3.
struct Interval {
    Bits rows, cols;
};

// ComputeIntervals: Algorithm 2 of Belohlavek & Trnecka (2015).
bool compute_intervals(const Context& I, const Context& E, const Clock& clk,
                       std::vector<Interval>& G) {
    Uncovered U(E.row);
    const Bits all_cols(I.n, true);
    const Bits full_rows = I.down(all_cols);  // C^up^down when (D u {j})^down_E is empty
    Bits a(I.m), q1(I.n), r(I.m), dd(I.n), p(I.m), q(I.n);
    Bits best_c, best_d, best_r, best_q;
    while (U.total > 0) {
        if (clk.over()) return false;
        Bits c(I.m, true), d(I.n);
        Bits rem_r(I.m), rem_q(I.n);
        long long v = 0;
        for (;;) {
            long long best = v;
            int best_j = -1;
            for (int j = 0; j < I.n; j++) {
                if (d.get(j)) continue;
                a.assign_and(c, E.col[j]);  // C = (D u {j})^down_E
                bool empty = a.none();
                if (empty) r = full_rows;
                else {
                    I.up(a, q1);
                    I.down(q1, r);          // C^up_I^down_I
                }
                if (U.row_sum(r) <= best) continue;
                if (empty) q = all_cols;
                else {
                    E.up(a, dd);            // D = C^up_E
                    I.down(dd, p);
                    I.up(p, q);             // D^down_I^up_I
                }
                long long val = U.area(r, q);
                if (val > best) {
                    best = val;
                    best_j = j;
                    best_c = a;
                    best_d = empty ? all_cols : dd;
                    best_r = r;
                    best_q = q;
                }
            }
            if (best_j < 0) break;
            c = best_c;
            d = best_d;
            v = best;
            rem_r = best_r;
            rem_q = best_q;
        }
        if (v == 0) break;  // unreachable: an uncovered essential entry always gives v > 0
        G.push_back({I.down(d), I.up(c)});
        U.remove(rem_r, rem_q);
    }
    return true;
}

// GreConD-like search inside the interval (lines 6-12 of Algorithm 1).
long long best_in_interval(const Context& I, const Uncovered& U, const Interval& g,
                           Bits& e, Bits& f) {
    e = g.rows;
    f = Bits(I.n);
    const long long ny = g.cols.count();
    long long s = 0;
    Bits a(I.m), b(I.n), best_e, best_f;
    for (;;) {
        long long best = s;
        int best_j = -1;
        for (int j = 0; j < I.n; j++) {
            if (!g.cols.get(j) || f.get(j)) continue;
            a.assign_and(e, I.col[j]);  // (F u {j})^down_J
            long long ca = a.count();
            if (ca == 0 || ca * ny <= best || U.row_sum(a) <= best) continue;
            I.up(a, b);
            b &= g.cols;                // ^up_J
            long long val = U.area(a, b);
            if (val > best) {
                best = val;
                best_j = j;
                best_e = a;
                best_f = b;
            }
        }
        if (best_j < 0) break;
        e = best_e;
        f = best_f;
        s = best;
    }
    return s;
}

// GreEss: Algorithm 1 of Belohlavek & Trnecka (2015).
void greess(const Context& I, State& st, int k, int& n_intervals, int& n_essential) {
    Context E = essential_part(I);
    n_essential = 0;
    for (const auto& r : E.row) n_essential += r.count();
    std::vector<Interval> G;
    bool done = compute_intervals(I, E, st.clk, G);
    n_intervals = (int)G.size();
    st.start();
    if (!done) { st.run.stop = "timeout_intervals"; return; }

    // The best factor of an interval depends only on U inside its rows x cols,
    // so it is recomputed only after a removal that touches that block.
    std::vector<char> used(G.size(), 0), dirty(G.size(), 1);
    std::vector<long long> val(G.size(), 0);
    std::vector<Bits> fe(G.size()), ff(G.size());
    size_t n_used = 0;
    bool exclusion = true;  // as in the reference code: lifted once every interval is used
    while (st.U.total > 0) {
        if ((int)st.run.factors.size() >= k) return;
        long long s = 0;
        int sel = -1;
        for (size_t g = 0; g < G.size(); g++) {
            if (exclusion && used[g]) continue;
            if (st.clk.over()) { st.run.stop = "timeout"; return; }
            if (dirty[g]) {
                val[g] = best_in_interval(I, st.U, G[g], fe[g], ff[g]);
                dirty[g] = 0;
            }
            if (val[g] > s) {
                s = val[g];
                sel = (int)g;
            }
        }
        if (sel < 0) { st.run.stop = "no_progress"; return; }
        Bits c = fe[sel], d = ff[sel];
        st.add(c, d);
        if (!used[sel]) {
            used[sel] = 1;
            if (++n_used == G.size()) exclusion = false;
        }
        for (size_t g = 0; g < G.size(); g++)
            if (c.intersects(G[g].rows) && d.intersects(G[g].cols)) dirty[g] = 1;
    }
    st.run.stop = "exact";
}

// -(c log2(c/t) + (t-c) log2((t-c)/t))
double code_bits(double c, double t) {
    double r = 0;
    if (c > 0) r -= c * std::log2(c / t);
    if (t - c > 0) r -= (t - c) * std::log2((t - c) / t);
    return r;
}

// MDLGreConD: Algorithm 1 of Makhalova & Trnecka (2021), cost (6) of Section 4.1.
void mdlgrecond(const Context& I, State& st, int k) {
    const double m = I.m, n = I.n;
    auto factor_bits = [&](long long ce, long long ci) {
        return std::log2(m) + code_bits(ce, m) + std::log2(n) + code_bits(ci, n);
    };
    auto total = [&](int nf, double fbits, long long err) {
        return 3 * std::log2(std::max({m, n, double(nf)})) + fbits + std::log2(m * n) +
               code_bits(err, m * n);
    };
    auto in_model = [&](const Bits& c) {
        for (const auto& f : st.run.factors) if (f.ext == c) return true;
        return false;
    };

    st.start();
    double total_cost = std::numeric_limits<double>::infinity();
    double fbits = 0;
    long long err = st.U.total;
    Bits a(I.m), b(I.n), best_c(I.m), best_d(I.n);
    for (;;) {
        if ((int)st.run.factors.size() >= k) return;
        if (st.clk.over()) { st.run.stop = "timeout"; return; }
        const int nf = (int)st.run.factors.size() + 1;
        Bits c(I.m, true), d(I.n);
        bool changed = false;
        long long gain = 0;
        for (;;) {
            double tc = total_cost;
            int best_j = -1;
            long long best_gain = 0;
            for (int j = 0; j < I.n; j++) {
                if (d.get(j)) continue;
                a.assign_and(c, I.col[j]);  // C' = D'^down with D' = (D u {j})^down^up
                I.up(a, b);
                if (in_model(a)) continue;
                long long g = st.U.area(a, b);
                double cost = total(nf, fbits + factor_bits(a.count(), b.count()), err - g);
                if (cost < tc) {
                    tc = cost;
                    best_j = j;
                    best_gain = g;
                    best_c = a;
                    best_d = b;
                }
            }
            total_cost = tc;
            if (best_j < 0) break;
            c = best_c;
            d = best_d;
            gain = best_gain;
            changed = true;
        }
        if (!changed) { st.run.stop = "mdl"; return; }
        fbits += factor_bits(c.count(), d.count());
        err -= gain;
        st.add(c, d);
    }
}

struct Eval {
    long long uncovered = 0, covered_zeros = 0;
};

// Independent check of A o B against the input matrix, on known cells.
Eval evaluate(const Matrix& M, const std::vector<Factor>& fs, size_t l) {
    Eval e;
    std::vector<Bits> prod(M.rows, Bits(M.cols));
    for (size_t f = 0; f < l; f++) fs[f].ext.for_each([&](int i) { prod[i] |= fs[f].itt; });
    for (int i = 0; i < M.rows; i++)
        for (int j = 0; j < M.cols; j++) {
            bool p = prod[i].get(j);
            if (M(i, j) == 1 && !p) e.uncovered++;
            if (M(i, j) == 0 && p) e.covered_zeros++;
        }
    return e;
}

void write_factors(const std::vector<Factor>& fs, size_t l, int m, int n,
                   const std::string& a_path, const std::string& b_path) {
    std::ofstream fa(a_path), fb(b_path);
    for (int i = 0; i < m; i++) {
        for (size_t f = 0; f < l; f++) fa << (f ? "," : "") << fs[f].ext.get(i);
        fa << "\n";
    }
    for (size_t f = 0; f < l; f++) {
        for (int j = 0; j < n; j++) fb << (j ? "," : "") << fs[f].itt.get(j);
        fb << "\n";
    }
}

void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <grecond|greess|mdlgrecond> <matrix.csv> <k>\n"
              << "       [--cuts k1,k2,...] [--missing free|zero|one] [--time-limit s] [--out-dir dir]\n"
              << "Factors are computed up to max(cuts) (default: k); one result row and one\n"
              << "A/B pair is written per cut, from the first cut factors.\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) { usage(argv[0]); return 2; }
    const std::string algo = argv[1], path = argv[2];
    const int k = std::stoi(argv[3]);
    std::vector<int> cuts;
    std::string missing = "free", out_dir = ".";
    Clock clk;
    for (int i = 4; i < argc; i++) {
        std::string opt = argv[i];
        if (i + 1 >= argc) { usage(argv[0]); return 2; }
        std::string val = argv[++i];
        if (opt == "--cuts") {
            std::stringstream ss(val);
            for (std::string t; std::getline(ss, t, ',');) cuts.push_back(std::stoi(t));
        } else if (opt == "--missing") missing = val;
        else if (opt == "--time-limit") clk.limit = std::stod(val);
        else if (opt == "--out-dir") out_dir = val;
        else { usage(argv[0]); return 2; }
    }
    if (cuts.empty()) cuts.push_back(k);
    if ((algo != "grecond" && algo != "greess" && algo != "mdlgrecond") ||
        (missing != "free" && missing != "zero" && missing != "one")) {
        usage(argv[0]);
        return 2;
    }

    Matrix M = CSVMatrixLoader::loadFromCSV(path);
    if (M.rows == 0) return 1;
    const int m = M.rows, n = M.cols;
    Context I(m, n);
    std::vector<Bits> known_ones(m, Bits(n));
    long long ones = 0, unknown = 0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            if (M(i, j) == 1) { known_ones[i].set(j); ones++; }
            if (M(i, j) == -1) unknown++;
            if (M(i, j) == 1 || (M(i, j) == -1 && missing != "zero")) I.set(i, j);
        }

    const int k_max = *std::max_element(cuts.begin(), cuts.end());
    State st(missing == "one" ? I.row : known_ones, known_ones, clk);
    int n_intervals = -1, n_essential = -1;
    clk.start = std::chrono::steady_clock::now();  // reading the CSV is not timed
    if (algo == "grecond") grecond(I, st, k_max);
    else if (algo == "greess") greess(I, st, k_max, n_intervals, n_essential);
    else mdlgrecond(I, st, k_max);
    const double run_time = clk.elapsed();
    const Run& run = st.run;

    std::filesystem::create_directories(out_dir);
    const std::string base = std::filesystem::path(path).stem().string();
    const std::string stem = out_dir + "/" + base + "_" + algo + "_" + missing;

    {
        std::ofstream mf(out_dir + "/" + base + "_M.csv");
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) mf << (j ? "," : "") << M(i, j);
            mf << "\n";
        }
        std::ofstream tf(stem + ".trace");
        tf << "factors,uncovered_ones,time_s\n";
        for (size_t l = 0; l < run.uncovered_known.size(); l++)
            tf << l << "," << run.uncovered_known[l] << "," << std::fixed << std::setprecision(6)
               << run.time_s[l] << "\n";
    }

    std::ofstream out(stem + ".csv");
    out << "dataset,algo,missing,k,factors,nonempty_factors,ones,unknown,uncovered_ones,covered_zeros,error,"
           "coverage,time_s,total_time_s,stop,essential_ones,intervals,self_check\n";
    bool ok = true;
    std::cout << "[" << algo << "] " << base << " " << m << "x" << n << " ones=" << ones
              << " unknown=" << unknown << " missing=" << missing << " factors="
              << run.factors.size() << " stop=" << run.stop << " time=" << run_time << "s\n";
    for (int cut : cuts) {
        size_t l = std::min<size_t>(cut, run.factors.size());
        if (run.uncovered_known.empty()) l = 0;
        Eval e = evaluate(M, run.factors, l);
        bool check = run.uncovered_known.size() > l && e.uncovered == run.uncovered_known[l] &&
                     e.covered_zeros == 0;
        ok = ok && check;
        double t = run.time_s.size() > l ? run.time_s[l] : run_time;
        std::string stop = l == run.factors.size() ? run.stop : "k_reached";
        size_t nonempty = 0;
        for (size_t f = 0; f < l; f++) nonempty += !run.factors[f].ext.none() && !run.factors[f].itt.none();
        out << base << "," << algo << "," << missing << "," << cut << "," << l << "," << nonempty << "," << ones << ","
            << unknown << "," << e.uncovered << "," << e.covered_zeros << ","
            << e.uncovered + e.covered_zeros << "," << std::setprecision(6)
            << (ones ? 1.0 - double(e.uncovered) / ones : 1.0) << "," << std::setprecision(6)
            << std::fixed << t << "," << run_time << std::defaultfloat << "," << stop << ","
            << n_essential << "," << n_intervals << "," << (check ? "ok" : "FAILED") << "\n";
        std::cout << "  k=" << cut << " factors=" << l << " error=" << e.uncovered + e.covered_zeros
                  << " covered_zeros=" << e.covered_zeros << " time=" << t << "s"
                  << (check ? "" : "  SELF-CHECK FAILED") << "\n";
        write_factors(run.factors, l, m, n, stem + "_k" + std::to_string(cut) + "_A.csv",
                      stem + "_k" + std::to_string(cut) + "_B.csv");
    }
    return ok ? 0 : 3;
}
