#pragma once
#include "Matrix.hpp"
#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <atomic>
#include <cmath>

// utilisee a la place du unordered_set 
class SparseIntSet {
    std::vector<int> items;
    std::vector<int> pos;
public:
    SparseIntSet() = default;
    void resize_capacity(int max_val) {
        if ((int)pos.size() < max_val) pos.resize(max_val, -1);
    }
    void clear() {
        for (int x : items) pos[x] = -1;
        items.clear();
    }
    void insert(int x) {
        if (pos[x] != -1) return;
        pos[x] = (int)items.size();
        items.push_back(x);
    }
    void erase(int x) {
        int p = pos[x];
        if (p == -1) return;
        int last = items.back();
        items[p] = last;
        pos[last] = p;
        items.pop_back();
        pos[x] = -1;
    }
    bool contains(int x) const { return x >= 0 && x < (int)pos.size() && pos[x] != -1; }
    int get_random(std::mt19937& rng) const {
        std::uniform_int_distribution<size_t> d(0, items.size() - 1);
        return items[d(rng)];
    }
    size_t size() const { return items.size(); }
    bool empty() const { return items.empty(); }
    const std::vector<int>& as_vector() const { return items; }
};

struct LocalSearchResult {
    bool success;
    double total_time;
    double last_improvement_ms;
    int iterations;
    int initial_errors;
    int final_errors;
    std::vector<int> error_history;
    Matrix A_solution;
    Matrix B_solution;
    std::string method_name;
    std::string stop_reason;

    LocalSearchResult() : success(false), total_time(0), last_improvement_ms(0),
                          iterations(0), initial_errors(0), final_errors(0),
                          A_solution(0, 0), B_solution(0, 0),
                          method_name("BASIC"), stop_reason("") {}
};

class BMFLocalSearch {
public:
    int m, n, k;
    Matrix M;
    Matrix A, B;

    std::vector<int> nbcover_flat;
    inline int& nbcover_ref(int i, int j)       { return nbcover_flat[i * n + j]; }
    inline int  nbcover(int i, int j) const     { return nbcover_flat[i * n + j]; }

    std::vector<SparseIntSet> rows_A; // rows_A[l] = {i : A[i,l]=1}
    std::vector<SparseIntSet> cols_B; // cols_B[l] = {j : B[l,j]=1}
    std::vector<SparseIntSet> cols_A; // cols_A[i] = {l : A[i,l]=1}
    std::vector<SparseIntSet> rows_B; // rows_B[j] = {l : B[l,j]=1}

    std::vector<SparseIntSet> zero_cover_in_row;
    std::vector<SparseIntSet> zero_cover_in_col;

    std::vector<double> score_A_flat;
    std::vector<double> score_B_flat;
    inline double& score_A_ref(int i, int l)    { return score_A_flat[i * k + l]; }
    inline double  score_A(int i, int l) const  { return score_A_flat[i * k + l]; }
    inline double& score_B_ref(int l, int j)    { return score_B_flat[l * n + j]; }
    inline double  score_B(int l, int j) const  { return score_B_flat[l * n + j]; }

    SparseIntSet pos_score_A;
    SparseIntSet pos_score_B;

    std::vector<double> weight_flat;
    inline double& weight_ref(int i, int j)     { return weight_flat[i * n + j]; }
    inline double  weight(int i, int j) const   { return weight_flat[i * n + j]; }

    SparseIntSet unsat_hard_cells; // M=0 et nbcover>=1
    SparseIntSet unsat_soft_cells; // M=1 et nbcover=0

    std::mt19937 rng;

    BMFLocalSearch(int k_factors, const Matrix& matrix, unsigned seed = 42);

    void initialize_greedy();
    void compute_all_counts();

    int count_errors() const;

    void flip_A(int i, int l);
    void flip_B(int l, int j);

    void apply_weight_change(int i, int j, double dw);

    LocalSearchResult solve_weighted(int max_iterations = 100000,
                                     bool verbose = false,
                                     std::atomic<bool>* stop_flag = nullptr,
                                     double timeout_ms = 0,
                                     double smooth_prob_override = -1.0,
                                     double h_inc = 1.0,
                                     double s_inc = 1.0);

    // M=1 (soft): erreur ssi nbcover=0. M=0 (hard): erreur ssi nbcover>=1. M=-1: ignore.
    static inline int err_val(int m_ij, int nb) {
        if (m_ij == -1) return 0;
        return m_ij ? (nb == 0 ? 1 : 0) : (nb >= 1 ? 1 : 0);
    }

    // Contribution au score (positif = ameliorant) d'une cellule pour un flip
    // qui changerait nb par sign in {+1, -1}.
    static inline double score_contrib(int m_ij, int nb, int sign, double w) {
        if (m_ij == -1) return 0.0;
        int de = err_val(m_ij, nb + sign) - err_val(m_ij, nb);
        return -w * (double)de;
    }

private:
    inline void update_pos_score_A(int i, int l) {
        int idx = i * k + l;
        bool good = score_A_flat[idx] > 1e-9;
        bool in   = pos_score_A.contains(idx);
        if (good && !in)  pos_score_A.insert(idx);
        if (!good && in)  pos_score_A.erase(idx);
    }
    inline void update_pos_score_B(int l, int j) {
        int idx = l * n + j;
        bool good = score_B_flat[idx] > 1e-9;
        bool in   = pos_score_B.contains(idx);
        if (good && !in)  pos_score_B.insert(idx);
        if (!good && in)  pos_score_B.erase(idx);
    }

    bool pick_good_flip(char& type, int& a, int& b);
    bool try_double_flip();
    void update_weights(double sp, double h_inc, double s_inc, bool& did_smooth);
    void force_flip_from_unsat();

    inline double effective_smooth_prob() const {
        long long nvars = (long long)m * k + (long long)k * n;
        return (nvars < 10000) ? 0.01 : 0.0003;
    }
};
