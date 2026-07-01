#pragma once
#include <vector>
#include <map>
#include <cstdint>

// Forward declaration - IPAMIR C API (defined in ipamir.h)
extern "C" {
    void* ipamir_init();
    void ipamir_release(void* solver);
    void ipamir_add_hard(void* solver, int32_t lit_or_zero);
    void ipamir_add_soft_lit(void* solver, int32_t lit, uint64_t weight);
    void ipamir_assume(void* solver, int32_t lit);
    int ipamir_solve(void* solver);
    uint64_t ipamir_val_obj(void* solver);
    int32_t ipamir_val_lit(void* solver, int32_t lit);
}

// IPAMIR return codes
constexpr int IPAMIR_RESULT_INTERRUPTED = 0;
constexpr int IPAMIR_RESULT_SAT = 10;
constexpr int IPAMIR_RESULT_UNSAT = 20;
constexpr int IPAMIR_RESULT_OPTIMAL = 30;
constexpr int IPAMIR_RESULT_ERROR = 40;

/**
 * @class SATSolver
 * @brief Wrapper for the EvalMaxSAT2022 MaxSAT solver via the IPAMIR API
 *
 * Provides a simplified interface for solving weighted MaxSAT problems.
 * The IPAMIR API natively supports assumptions, incremental clause addition,
 * and weighted soft clauses.
 */
class SATSolver {
private:
    void* ipamir_solver;

public:
    int num_vars;
    std::vector<bool> model;

    SATSolver();
    ~SATSolver();

    /** @brief Creates a new SAT variable. @return The variable identifier (1-indexed). */
    int new_var();

    /** @brief Adds a hard clause. @param clause Vector of literals. */
    void add_clause(const std::vector<int>& clause);

    /**
     * @brief Adds a weighted soft clause.
     * @param clause Vector of literals.
     * @param weight Clause weight.
     * @return The relaxation variable identifier (for non-unit clauses).
     *
     * For unit clauses, uses ipamir_add_soft_lit directly.
     * For non-unit clauses, introduces an auxiliary variable.
     */
    int add_soft_clause(const std::vector<int>& clause, long long weight = 1);

    /** @brief Solves the MaxSAT problem. @return true if a solution was found (SAT or OPTIMAL). */
    bool solve();

    /**
     * @brief Solves the MaxSAT problem with native assumptions.
     * @param assumptions Vector of literals to assume (positive = true, negative = false).
     * @return true if a solution was found.
     *
     * Assumptions are automatically cleared after solve().
     */
    bool solve_with_assumptions(const std::vector<int>& assumptions);

    /** @brief Adds an assumption for the next solve() call. Cleared after each solve(). */
    void assume(int lit);

    /**
     * @brief Adds a conditional control clause (for incremental LNS).
     * Adds: ctrl_var -> target_lit (equivalent to: -ctrl_var | target_lit)
     */
    void add_control_clause(int ctrl_var, int target_lit);

    /** @brief Gets the value of a literal in the current model. */
    bool getValue(int lit);

    /** @brief Gets the solution cost (sum of weights of unsatisfied soft clauses). */
    long long get_cost();

    /** @brief Resets the solver (release + init). */
    void reset();

    /**
     * @brief Adds a conditional clause (with selector).
     * Adds: selector -> (clause), i.e., -selector | l1 | l2 | ... | ln
     */
    void add_conditional_clause(int selector, const std::vector<int>& clause);

    /**
     * @brief Adds a conditional unit constraint.
     * Adds: selector -> lit, i.e., -selector | lit
     */
    void add_conditional_unit(int selector, int lit);
};

/**
 * @class VariableManager
 * @brief Manages SAT variables indexed by matrix positions (i, j).
 */
class VariableManager {
public:
    std::map<std::pair<int,int>, int> vars;
    SATSolver& solver;

    VariableManager(SATSolver& s) : solver(s) {}

    /** @brief Gets or creates a variable for position (i, j). */
    int get(int i, int j) {
        auto key = std::make_pair(i, j);
        if (vars.find(key) == vars.end()) {
            vars[key] = solver.new_var();
        }
        return vars[key];
    }

    void reset() {
        vars.clear();
    }
};
