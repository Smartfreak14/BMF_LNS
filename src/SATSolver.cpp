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
        // IPAMIR: ajouter chaque littéral puis terminer par 0
        for (int lit : clause) {
            ipamir_add_hard(ipamir_solver, lit);
        }
        ipamir_add_hard(ipamir_solver, 0);  // Terminer la clause
    }
    clauses.push_back(clause);
}

int SATSolver::add_soft_clause(const std::vector<int>& clause, long long weight) {
    if (!ipamir_solver || clause.empty()) {
        return -1;
    }
    
    if (clause.size() == 1) {
        // Clause unitaire: utiliser directement ipamir_add_soft_lit
        // soft(l) avec poids w signifie: si l est VRAI, on paie w
        // Donc pour encourager la satisfaction de la clause {l}, 
        // on déclare -l comme soft literal (payer si -l est vrai = l est faux)
        ipamir_add_soft_lit(ipamir_solver, -clause[0], static_cast<uint64_t>(weight));
        return clause[0];
    } else {
        // Clause non-unitaire: introduire une variable auxiliaire b
        // Ajouter (C ∨ b) comme clause hard
        // Déclarer b comme soft literal avec poids weight
        // Si C n'est pas satisfaite, b doit être vrai → coût = weight
        int aux_var = new_var();
        
        // Ajouter (clause ∨ aux_var) comme hard
        for (int lit : clause) {
            ipamir_add_hard(ipamir_solver, lit);
        }
        ipamir_add_hard(ipamir_solver, aux_var);
        ipamir_add_hard(ipamir_solver, 0);
        
        // Déclarer aux_var comme soft (si aux_var=true, on paie weight)
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
        // Extraire le modèle
        model.clear();
        model.resize(num_vars + 1, false);
        for (int i = 1; i <= num_vars; i++) {
            int val = ipamir_val_lit(ipamir_solver, i);
            model[i] = (val > 0);  // val > 0 signifie que le littéral est vrai
        }
        return true;
    }
    
    return false;
}

bool SATSolver::solve_with_assumptions(const std::vector<int>& assumptions) {
    if (!ipamir_solver) {
        return false;
    }
    
    // Ajouter les assumptions via l'API IPAMIR native
    for (int lit : assumptions) {
        ipamir_assume(ipamir_solver, lit);
    }
    
    int result = ipamir_solve(ipamir_solver);
    // Note: les assumptions sont automatiquement effacées après solve()
    
    if (result == IPAMIR_RESULT_OPTIMAL || result == IPAMIR_RESULT_SAT) {
        // Extraire le modèle
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
        // ctrl_var → target_lit  ≡  ¬ctrl_var ∨ target_lit
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

void SATSolver::setIncremental(bool value) {
    // IPAMIR est toujours incrémental par conception
    // Cette fonction est conservée pour la compatibilité API
    (void)value;
}

void SATSolver::reset() {
    if (ipamir_solver) {
        ipamir_release(ipamir_solver);
    }
    ipamir_solver = ipamir_init();
    num_vars = 0;
    clauses.clear();
    model.clear();
}

void SATSolver::add_conditional_clause(int selector, const std::vector<int>& clause) {
    if (ipamir_solver && !clause.empty()) {
        // selector → (l1 ∨ l2 ∨ ... ∨ ln)  ≡  ¬selector ∨ l1 ∨ l2 ∨ ... ∨ ln
        ipamir_add_hard(ipamir_solver, -selector);
        for (int lit : clause) {
            ipamir_add_hard(ipamir_solver, lit);
        }
        ipamir_add_hard(ipamir_solver, 0);
    }
}

void SATSolver::add_conditional_unit(int selector, int lit) {
    if (ipamir_solver) {
        // selector → lit  ≡  ¬selector ∨ lit
        ipamir_add_hard(ipamir_solver, -selector);
        ipamir_add_hard(ipamir_solver, lit);
        ipamir_add_hard(ipamir_solver, 0);
    }
}
