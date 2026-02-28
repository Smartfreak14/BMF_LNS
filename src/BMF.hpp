#pragma once
#include "Matrix.hpp"
#include "SATSolver.hpp"
#include <chrono>
#include <functional>
#include <map>
#include <vector>
#include <climits>
#include <string>
#include <fstream>
#include <atomic>

// Résultat de factorisation BMF
struct BMFResult {
    bool success;           // true si M = A × B exactement
    double total_time;      // Temps en ms
    int k;                  // Nombre de facteurs utilisés
    int errors;             // Nombre d'erreurs d'alignement (0 si succès)
    Matrix A_solution;
    Matrix B_solution;
    
    BMFResult() : success(false), total_time(0), k(0), errors(0),
                  A_solution(0, 0), B_solution(0, 0) {}
};

// Résultat LNS pour BMF
struct LNSResult {
    bool success;
    double total_time;
    int iterations;
    int initial_errors;
    int final_errors;
    std::vector<int> error_history;
    std::vector<int> modification_history;
    std::vector<double> iteration_times;
    Matrix A_solution;
    Matrix B_solution;
    std::string init_method;
    double solver_init_time;
    int solver_resets;
    int total_assumptions;
    int stagnation_count;
    int max_stagnation;
    int total_stagnations;
    std::string stop_reason;
    
    LNSResult() : success(false), total_time(0), iterations(0),
                  initial_errors(0), final_errors(0),
                  A_solution(0, 0), B_solution(0, 0), init_method("random"),
                  solver_init_time(0), solver_resets(0), total_assumptions(0),
                  stagnation_count(0), max_stagnation(0), total_stagnations(0),
                  stop_reason("max_iterations") {}
};

// Phase d'optimisation pour LNS-v3
enum class OptimizationPhase {
    FIX_A,
    FIX_B,
    JOINT
};

class BMF {
public:
    int k, m, n;
    Matrix M;
    SATSolver solver;
    VariableManager A, B;
    
    BMF(int k_factors, const Matrix& matrix) 
        : k(k_factors), M(matrix), A(solver), B(solver),
          incremental_initialized(false), 
          initial_A(0, 0), initial_B(0, 0), has_initial_solution(false) {
        m = M.rows;
        n = M.cols;
    }
    
    std::pair<Matrix, Matrix> generate_random_AB(double density = 0.3, unsigned seed = 42);
    
    void set_initial_solution(const Matrix& A_init, const Matrix& B_init) {
        initial_A = A_init;
        initial_B = B_init;
        has_initial_solution = true;
    }
    
    int count_alignment_errors(const Matrix& A_mat, const Matrix& B_mat);
    std::pair<Matrix, Matrix> get_solution();
    std::pair<std::vector<int>, std::vector<int>> compute_errors(
        const Matrix& A_mat, const Matrix& B_mat);
    std::vector<int> select_top_k(const std::vector<int>& errors, int k);
    
    // Méthode LNS-v3
    LNSResult solve_lns_v3(
        int max_iterations = 200,
        double neighborhood_size = 0.2,
        unsigned seed = 42,
        int stagnation_threshold = 15,
        bool verbose = true,
        double timeout_ms = 0,
        std::atomic<bool>* stop_flag = nullptr
    );
    
    // Méthodes LNS publiques pour utilisation directe
    int lns_step_partial(
        Matrix& A_current, Matrix& B_current,
        const std::vector<int>& rows_to_free,
        const std::vector<int>& cols_to_free
    );
    
    int lns_step_fix_A(
        Matrix& A_current, Matrix& B_current,
        const std::vector<int>& cols_to_free
    );
    
    int lns_step_fix_B(
        Matrix& A_current, Matrix& B_current,
        const std::vector<int>& rows_to_free
    );
    
private:
    bool incremental_initialized;
    Matrix initial_A;
    Matrix initial_B;
    bool has_initial_solution;
    
    void reset() {
        solver.reset();
        A.reset();
        B.reset();
        incremental_initialized = false;
    }

    std::vector<std::vector<int>> is_one(int i, int j);
    std::vector<std::vector<int>> is_zero(int i, int j);
    std::vector<std::vector<int>> factorization();
    
    void add_hard_clauses(const std::vector<std::vector<int>>& clauses);
    void add_soft_clauses(const std::vector<std::vector<int>>& clauses, int weight);
    
    void initialize_incremental_solver();
    
    int lns_step_ultra_targeted(
        Matrix& A_current, Matrix& B_current
    );
    
    void extract_solution(Matrix& A_out, Matrix& B_out);
};