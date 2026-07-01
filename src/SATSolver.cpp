#include "SATSolver.hpp"
#include <iostream>
#include <algorithm>

SATSolver::SATSolver() : ipamir_solver(nullptr), num_vars(0) {
    ipamir_solver = ipamir_init();
}

SATSolver::~SATSolver() {
    if (ipamir_solver) {
        ipamir_release(ipamir_solver);
        ipamir_solver = nullptr;
    }
}

int SATSolver::new_var() {
    num_vars++;
    return num_vars;
}

void SATSolver::add_clause(const std::vector<int>& clause) {
    if (ipamir_solver && !clause.empty()) {
        for (int lit : clause) {
            ipamir_add_hard(ipamir_solver, lit);
        }
        ipamir_add_hard(ipamir_solver, 0);
    }
}

int SATSolver::add_soft_clause(const std::vector<int>& clause, long long weight) {
    if (!ipamir_solver || clause.empty()) {
        return -1;
    }

    if (clause.size() == 1) {
        ipamir_add_soft_lit(ipamir_solver, -clause[0], static_cast<uint64_t>(weight));
        return clause[0];
    } else {
        int aux_var = new_var();

        for (int lit : clause) {
            ipamir_add_hard(ipamir_solver, lit);
        }
        ipamir_add_hard(ipamir_solver, aux_var);
        ipamir_add_hard(ipamir_solver, 0);

        ipamir_add_soft_lit(ipamir_solver, aux_var, static_cast<uint64_t>(weight));
        return aux_var;
    }
}

bool SATSolver::solve() {
    if (!ipamir_solver) {
        return false;
    }

    int result = ipamir_solve(ipamir_solver);

    if (result == IPAMIR_RESULT_OPTIMAL || result == IPAMIR_RESULT_SAT) {
        model.clear();
        model.resize(num_vars + 1, false);
        for (int i = 1; i <= num_vars; i++) {
            int val = ipamir_val_lit(ipamir_solver, i);
            model[i] = (val > 0);
        }
        return true;
    }

    return false;
}

bool SATSolver::solve_with_assumptions(const std::vector<int>& assumptions) {
    if (!ipamir_solver) {
        return false;
    }

    for (int lit : assumptions) {
        ipamir_assume(ipamir_solver, lit);
    }

    int result = ipamir_solve(ipamir_solver);

    if (result == IPAMIR_RESULT_OPTIMAL || result == IPAMIR_RESULT_SAT) {
        model.clear();
        model.resize(num_vars + 1, false);
        for (int i = 1; i <= num_vars; i++) {
            int val = ipamir_val_lit(ipamir_solver, i);
            model[i] = (val > 0);
        }
        return true;
    }

    return false;
}

void SATSolver::assume(int lit) {
    if (ipamir_solver) {
        ipamir_assume(ipamir_solver, lit);
    }
}

void SATSolver::add_control_clause(int ctrl_var, int target_lit) {
    if (ipamir_solver) {
        ipamir_add_hard(ipamir_solver, -ctrl_var);
        ipamir_add_hard(ipamir_solver, target_lit);
        ipamir_add_hard(ipamir_solver, 0);
    }
}

bool SATSolver::getValue(int lit) {
    if (!ipamir_solver) {
        return false;
    }
    int val = ipamir_val_lit(ipamir_solver, lit);
    return (val > 0);
}

long long SATSolver::get_cost() {
    if (!ipamir_solver) {
        return -1;
    }
    return static_cast<long long>(ipamir_val_obj(ipamir_solver));
}

void SATSolver::reset() {
    if (ipamir_solver) {
        ipamir_release(ipamir_solver);
    }
    ipamir_solver = ipamir_init();
    num_vars = 0;
    model.clear();
}

void SATSolver::add_conditional_clause(int selector, const std::vector<int>& clause) {
    if (ipamir_solver && !clause.empty()) {
        ipamir_add_hard(ipamir_solver, -selector);
        for (int lit : clause) {
            ipamir_add_hard(ipamir_solver, lit);
        }
        ipamir_add_hard(ipamir_solver, 0);
    }
}

void SATSolver::add_conditional_unit(int selector, int lit) {
    if (ipamir_solver) {
        ipamir_add_hard(ipamir_solver, -selector);
        ipamir_add_hard(ipamir_solver, lit);
        ipamir_add_hard(ipamir_solver, 0);
    }
}
