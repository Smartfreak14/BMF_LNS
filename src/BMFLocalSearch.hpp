#pragma once
#include "Matrix.hpp"
#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <unordered_set>
#include <atomic>

struct LocalSearchResult {
    bool success;
    double total_time;
    int iterations;
    int initial_errors;
    int final_errors;
    std::vector<int> error_history;
    Matrix A_solution;
    Matrix B_solution;
    std::string method_name;
    std::string stop_reason;
    
    LocalSearchResult() : success(false), total_time(0), iterations(0),
                          initial_errors(0), final_errors(0),
                          A_solution(0, 0), B_solution(0, 0),
                          method_name("BASIC"), stop_reason("") {}
};

class BMFLocalSearch {
public:
    int m, n, k;
    Matrix M;
    Matrix A, B;
    
    std::vector<int> count_flat;
    std::vector<int> score_A_flat;
    std::vector<int> score_B_flat;
    
    inline int& count_at(int i, int j) { return count_flat[i * n + j]; }
    inline int count_at(int i, int j) const { return count_flat[i * n + j]; }
    inline int& score_A_at(int i, int l) { return score_A_flat[i * k + l]; }
    inline int score_A_at(int i, int l) const { return score_A_flat[i * k + l]; }
    inline int& score_B_at(int l, int j) { return score_B_flat[l * n + j]; }
    inline int score_B_at(int l, int j) const { return score_B_flat[l * n + j]; }
    
    std::vector<std::unordered_set<int>> A_ones_set;
    std::vector<std::unordered_set<int>> B_ones_set;
    
    std::vector<std::vector<int>> count;
    std::vector<std::vector<int>> score_A;
    std::vector<std::vector<int>> score_B;
    std::vector<std::vector<int>> B_ones_cols;
    std::vector<std::vector<int>> A_ones_rows;
    
    std::mt19937 rng;
    
    BMFLocalSearch(int k_factors, const Matrix& matrix, unsigned seed = 42);
    
    void initialize_greedy();
    void compute_all_counts();
    void flip_A(int i, int l);
    void flip_B(int l, int j);
    int count_errors() const;
    
    void compute_all_scores();
    int compute_score_A(int i, int l);
    int compute_score_B(int l, int j);
    std::tuple<char, int, int, int> find_best_flip();
    LocalSearchResult solve_basic(int max_iterations = 100000, bool verbose = false);

    LocalSearchResult solve_weighted(int max_iterations = 100000,
                                     double penalty_increment = 0.1,
                                     int max_stagnation = 100,
                                     bool verbose = false,
                                     std::atomic<bool>* stop_flag = nullptr);
    LocalSearchResult solve(int max_iterations = 100000, bool verbose = false);
    
    std::vector<double> weights_flat;
    inline double& weight_at(int i, int j) { return weights_flat[i * n + j]; }
    inline double weight_at(int i, int j) const { return weights_flat[i * n + j]; }
    
    double compute_weighted_score_A(int i, int l);
    double compute_weighted_score_B(int l, int j);
    std::tuple<char, int, int, double> find_best_weighted_flip(int current_errors);
};
