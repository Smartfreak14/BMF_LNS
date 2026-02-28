#pragma once
#include <vector>
#include <map>
#include <cstdint>

// Forward declaration - IPAMIR C API (défini dans ipamir.h)
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

// Codes de retour IPAMIR
constexpr int IPAMIR_RESULT_INTERRUPTED = 0;
constexpr int IPAMIR_RESULT_SAT = 10;
constexpr int IPAMIR_RESULT_UNSAT = 20;
constexpr int IPAMIR_RESULT_OPTIMAL = 30;
constexpr int IPAMIR_RESULT_ERROR = 40;

/**
 * @class SATSolver
 * @brief Wrapper pour le solveur MaxSAT EvalMaxSAT2022 via l'API IPAMIR
 * 
 * Cette classe encapsule EvalMaxSAT2022 via l'interface IPAMIR pour fournir
 * une interface simplifiée pour résoudre des problèmes MaxSAT pondérés.
 * 
 * L'API IPAMIR supporte nativement:
 * - Les assumptions (ipamir_assume)
 * - L'ajout incrémental de clauses
 * - Les clauses soft pondérées
 */
class SATSolver {
private:
    void* ipamir_solver;

public:
    int num_vars;
    std::vector<std::vector<int>> clauses;  // Pour debug/tracking
    std::vector<bool> model;
    
    SATSolver();
    ~SATSolver();
    
    /**
     * @brief Crée une nouvelle variable SAT
     * @return L'identifiant de la nouvelle variable (1-indexed)
     */
    int new_var();
    
    /**
     * @brief Ajoute une clause dure (hard clause)
     * @param clause Vecteur de littéraux
     */
    void add_clause(const std::vector<int>& clause);
    
    /**
     * @brief Ajoute une clause souple (soft clause) avec un poids
     * @param clause Vecteur de littéraux
     * @param weight Poids de la clause
     * @return L'identifiant de la variable de relaxation (pour clauses non-unitaires)
     * 
     * Note: Pour les clauses unitaires, utilise directement ipamir_add_soft_lit.
     * Pour les clauses non-unitaires, introduit une variable auxiliaire.
     */
    int add_soft_clause(const std::vector<int>& clause, long long weight = 1);
    
    /**
     * @brief Résout le problème MaxSAT
     * @return true si une solution a été trouvée (SAT ou OPTIMAL)
     */
    bool solve();
    
    /**
     * @brief Résout le problème MaxSAT avec des assumptions natives
     * @param assumptions Vecteur de littéraux à assumer (positif = vrai, négatif = faux)
     * @return true si une solution a été trouvée
     * 
     * Note: Les assumptions sont automatiquement effacées après solve().
     * C'est la vraie implémentation incrémentale via IPAMIR.
     */
    bool solve_with_assumptions(const std::vector<int>& assumptions);
    
    /**
     * @brief Ajoute une assumption pour le prochain appel à solve()
     * @param lit Le littéral à assumer (positif = vrai, négatif = faux)
     * 
     * Les assumptions sont effacées après chaque appel à solve().
     */
    void assume(int lit);
    
    /**
     * @brief Ajoute une clause de contrôle conditionnelle (pour LNS incrémental)
     * @param ctrl_var Variable de contrôle (si vraie, la clause target doit être satisfaite)
     * @param target_lit Littéral cible à forcer quand ctrl_var est vrai
     * 
     * Ajoute: ctrl_var → target_lit (équivalent à: ¬ctrl_var ∨ target_lit)
     */
    void add_control_clause(int ctrl_var, int target_lit);
    
    /**
     * @brief Obtient la valeur d'un littéral dans le modèle
     * @param lit Le littéral (positif ou négatif)
     * @return La valeur booléenne du littéral
     */
    bool getValue(int lit);
    
    /**
     * @brief Obtient le coût de la solution (somme des poids des soft clauses non satisfaites)
     * @return Le coût
     */
    long long get_cost();
    
    /**
     * @brief Active le mode incrémental (no-op avec IPAMIR, toujours incrémental)
     * @param value Ignoré - IPAMIR est toujours incrémental
     */
    void setIncremental(bool value = true);
    
    /**
     * @brief Réinitialise le solveur (release + init)
     */
    void reset();
    
    /**
     * @brief Ajoute une clause conditionnelle (avec sélecteur)
     * @param selector Variable de sélection (si vraie, la clause doit être satisfaite)
     * @param clause Clause à ajouter conditionnellement
     * 
     * Ajoute: selector → (clause)  ≡  ¬selector ∨ clause
     */
    void add_conditional_clause(int selector, const std::vector<int>& clause);
    
    /**
     * @brief Ajoute une contrainte unitaire conditionnelle
     * @param selector Variable de sélection
     * @param lit Littéral à forcer quand selector est vrai
     * 
     * Ajoute: selector → lit  ≡  ¬selector ∨ lit
     */
    void add_conditional_unit(int selector, int lit);
};

/**
 * @class VariableManager
 * @brief Gestionnaire de variables SAT pour matrices indexées par (i, j)
 * 
 * Permet de créer et récupérer des variables SAT associées à des positions matricielles.
 */
class VariableManager {
public:
    std::map<std::pair<int,int>, int> vars;
    SATSolver& solver;
    
    VariableManager(SATSolver& s) : solver(s) {}
    
    /**
     * @brief Obtient ou crée une variable pour la position (i, j)
     * @param i Index de ligne
     * @param j Index de colonne
     * @return L'identifiant de la variable
     */
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
