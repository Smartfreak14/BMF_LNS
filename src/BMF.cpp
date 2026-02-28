#include "BMF.hpp"
#include "BMFLocalSearch.hpp"
#include <iostream>
#include <algorithm>
#include <random>
#include <numeric>
#include <cmath>
#include <limits>
#include <set>
#include <fstream>

// ==================== Contraintes de base ====================

std::vector<std::vector<int>> BMF::is_one(int i, int j) {
    std::vector<std::vector<int>> clauses;
    std::vector<int> disjunction;
    
    for (int l = 0; l < k; l++) {
        int t_var = solver.new_var();
        disjunction.push_back(t_var);
        
        clauses.push_back({-t_var, A.get(i, l)});
        clauses.push_back({-t_var, B.get(l, j)});
        clauses.push_back({-A.get(i, l), -B.get(l, j), t_var});
    }
    
    clauses.push_back(disjunction);
    return clauses;
}

std::vector<std::vector<int>> BMF::is_zero(int i, int j) {
    std::vector<std::vector<int>> clauses;
    
    for (int l = 0; l < k; l++) {
        clauses.push_back({-A.get(i, l), -B.get(l, j)});
    }
    
    return clauses;
}

std::vector<std::vector<int>> BMF::factorization() {
    std::vector<std::vector<int>> clauses;
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (M(i, j) == 1) {
                auto is_one_clauses = is_one(i, j);
                clauses.insert(clauses.end(), is_one_clauses.begin(), is_one_clauses.end());
            } else if (M(i, j) == 0) {
                auto is_zero_clauses = is_zero(i, j);
                clauses.insert(clauses.end(), is_zero_clauses.begin(), is_zero_clauses.end());
            }
        }
    }
    
    return clauses;
}

void BMF::add_hard_clauses(const std::vector<std::vector<int>>& clauses) {
    for (const auto& clause : clauses) {
        solver.add_clause(clause);
    }
}

void BMF::add_soft_clauses(const std::vector<std::vector<int>>& clauses, int weight) {
    for (const auto& clause : clauses) {
        solver.add_soft_clause(clause, weight);
    }
}

// ==================== Utilitaires ====================

std::pair<Matrix, Matrix> BMF::generate_random_AB(double density, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    Matrix A_rand(m, k, 0);
    Matrix B_rand(k, n, 0);
    
    for (int i = 0; i < m; i++) {
        for (int l = 0; l < k; l++) {
            if (dis(rng) < density) {
                A_rand(i, l) = 1;
            }
        }
    }
    
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            if (dis(rng) < density) {
                B_rand(l, j) = 1;
            }
        }
    }
    
    return {A_rand, B_rand};
}

int BMF::count_alignment_errors(const Matrix& A_mat, const Matrix& B_mat) {
    Matrix computed = A_mat.multiply(B_mat);
    int errors = 0;
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (M(i, j) == -1) continue;
            if (computed(i, j) != M(i, j)) {
                errors++;
            }
        }
    }
    
    return errors;
}

// ==================== SÉLECTION GUIDÉE PAR L'ERREUR ====================

std::pair<std::vector<int>, std::vector<int>> BMF::compute_errors(
    const Matrix& A_mat, const Matrix& B_mat
) {
    std::vector<int> row_errors(m, 0);
    std::vector<int> col_errors(n, 0);
    
    // Calcul direct SANS multiplication matricielle complète
    // Pour chaque cellule (i,j), calculer (A⊗B)[i,j] à la volée
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (M(i, j) == -1) continue;  // Don't care
            
            // Calcul (A⊗B)[i,j] = OR sur l : A[i,l] AND B[l,j]
            int product_ij = 0;
            for (int l = 0; l < k; l++) {
                if (A_mat(i, l) == 1 && B_mat(l, j) == 1) {
                    product_ij = 1;
                    break;  // OR booléen: dès qu'on trouve un 1, on arrête
                }
            }
            
            if (product_ij != M(i, j)) {
                row_errors[i]++;
                col_errors[j]++;
            }
        }
    }
    
    return {row_errors, col_errors};
}

std::vector<int> BMF::select_top_k(const std::vector<int>& errors, int k_select) {
    int n_items = errors.size();
    k_select = std::min(k_select, n_items);
    
    // Créer un vecteur d'indices triés par erreur décroissante
    std::vector<int> indices(n_items);
    std::iota(indices.begin(), indices.end(), 0);
    
    // Tri partiel: on ne trie que les k premiers éléments
    std::partial_sort(indices.begin(), indices.begin() + k_select, indices.end(),
        [&errors](int a, int b) {
            return errors[a] > errors[b];  // Tri décroissant
        });
    
    // Retourner les k premiers
    return std::vector<int>(indices.begin(), indices.begin() + k_select);
}

std::pair<Matrix, Matrix> BMF::get_solution() {
    if (solver.model.empty()) {
        return {Matrix(m, k, 0), Matrix(k, n, 0)};
    }
    
    Matrix A_sol(m, k, 0);
    for (int i = 0; i < m; i++) {
        for (int l = 0; l < k; l++) {
            auto it = A.vars.find({i, l});
            if (it != A.vars.end() && solver.getValue(it->second)) {
                A_sol(i, l) = 1;
            }
        }
    }
    
    Matrix B_sol(k, n, 0);
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            auto it = B.vars.find({l, j});
            if (it != B.vars.end() && solver.getValue(it->second)) {
                B_sol(l, j) = 1;
            }
        }
    }
    
    return {A_sol, B_sol};
}

void BMF::extract_solution(Matrix& A_out, Matrix& B_out) {
    for (int i = 0; i < m; i++) {
        for (int l = 0; l < k; l++) {
            auto it = A.vars.find({i, l});
            if (it != A.vars.end()) {
                A_out(i, l) = solver.getValue(it->second) ? 1 : 0;
            }
        }
    }
    
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            auto it = B.vars.find({l, j});
            if (it != B.vars.end()) {
                B_out(l, j) = solver.getValue(it->second) ? 1 : 0;
            }
        }
    }
}



// ==================== MODE INCRÉMENTAL  ====================

void BMF::initialize_incremental_solver() {
    if (incremental_initialized) return;
    
    // Créer toutes les variables A et B
    for (int i = 0; i < m; i++) {
        for (int l = 0; l < k; l++) {
            A.get(i, l);
        }
    }
    
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            B.get(l, j);
        }
    }
    
    // Ajouter toutes les clauses de factorisation comme SOFT
    // Ces clauses restent dans le solveur pour toutes les itérations
    add_soft_clauses(factorization(), 1);
    
    incremental_initialized = true;
}

// ==================== MODE ULTRA-CIBLÉ (pour peu d'erreurs) ====================

int BMF::lns_step_ultra_targeted(
    Matrix& A_current, Matrix& B_current
) {
    // Cette fonction est optimisée pour quand il reste TRÈS PEU d'erreurs (< 15)
    // Version améliorée: libère les lignes/colonnes d'erreur + leurs voisins factoriels
    
    // 1. Identifier les erreurs et les lignes/colonnes impliquées
    Matrix computed = A_current.multiply(B_current);
    std::set<int> error_rows_set, error_cols_set;
    std::vector<std::pair<int, int>> errors;  // Liste des cellules en erreur
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (M(i, j) != -1 && M(i, j) != computed(i, j)) {
                error_rows_set.insert(i);
                error_cols_set.insert(j);
                errors.push_back({i, j});
            }
        }
    }
    
    if (errors.empty()) return 0;  // Pas d'erreur !
    
    // 2. EXTENSION: Ajouter les lignes/colonnes qui partagent des facteurs actifs avec les erreurs
    std::set<int> extended_rows = error_rows_set;
    std::set<int> extended_cols = error_cols_set;
    
    // Trouver les facteurs impliqués dans les erreurs
    std::set<int> error_factors;
    for (auto& [i, j] : errors) {
        for (int l = 0; l < k; l++) {
            if (A_current(i, l) == 1 || B_current(l, j) == 1) {
                error_factors.insert(l);
            }
        }
    }
    
    // Ajouter les lignes/colonnes qui utilisent ces facteurs (limité pour performance)
    // AUGMENTÉ pour couvrir plus de contexte et sortir des minima locaux
    const int MAX_EXTRA_ROWS = 80;  // Maximum de lignes supplémentaires (augmenté de 50->80)
    const int MAX_EXTRA_COLS = 80;  // Maximum de colonnes supplémentaires (augmenté de 50->80)
    
    std::vector<std::pair<int, int>> row_factor_counts;  // (count, row)
    std::vector<std::pair<int, int>> col_factor_counts;  // (count, col)
    
    for (int i = 0; i < m; i++) {
        if (error_rows_set.count(i)) continue;  // Déjà incluse
        int count = 0;
        for (int l : error_factors) {
            if (A_current(i, l) == 1) count++;
        }
        if (count > 0) row_factor_counts.push_back({count, i});
    }
    
    for (int j = 0; j < n; j++) {
        if (error_cols_set.count(j)) continue;  // Déjà incluse
        int count = 0;
        for (int l : error_factors) {
            if (B_current(l, j) == 1) count++;
        }
        if (count > 0) col_factor_counts.push_back({count, j});
    }
    
    // Trier par nombre de facteurs partagés (décroissant)
    std::sort(row_factor_counts.begin(), row_factor_counts.end(), std::greater<>());
    std::sort(col_factor_counts.begin(), col_factor_counts.end(), std::greater<>());
    
    // Ajouter les top lignes/colonnes
    for (int i = 0; i < std::min((int)row_factor_counts.size(), MAX_EXTRA_ROWS); i++) {
        extended_rows.insert(row_factor_counts[i].second);
    }
    for (int j = 0; j < std::min((int)col_factor_counts.size(), MAX_EXTRA_COLS); j++) {
        extended_cols.insert(col_factor_counts[j].second);
    }
    
    std::vector<int> error_rows(extended_rows.begin(), extended_rows.end());
    std::vector<int> error_cols(extended_cols.begin(), extended_cols.end());
    
    std::cout << "  [ULTRA-CIBLÉ] " << errors.size() << " erreurs, "
              << error_rows.size() << " lignes × " << error_cols.size() << " colonnes"
              << " -> " << (error_rows.size() * k + k * error_cols.size()) << " variables"
              << " (+" << (error_rows.size() - error_rows_set.size()) << " lignes voisines, "
              << "+" << (error_cols.size() - error_cols_set.size()) << " colonnes voisines)" << std::endl;
    
    // 2. Créer le solveur avec uniquement les variables nécessaires
    SATSolver ultra_solver;
    VariableManager A_vars(ultra_solver), B_vars(ultra_solver);
    
    // Variables libres: A[error_rows, :] (toutes les colonnes de A pour ces lignes)
    for (int i : error_rows) {
        for (int l = 0; l < k; l++) {
            A_vars.get(i, l);
        }
    }
    
    // Variables libres: B[:, error_cols] (toutes les lignes de B pour ces colonnes)
    for (int l = 0; l < k; l++) {
        for (int j : error_cols) {
            B_vars.get(l, j);
        }
    }
    
    int hard_count = 0, soft_count = 0;
    
    // 3. Pour chaque cellule POTENTIELLEMENT affectée, créer les contraintes
    // Une cellule (i,j) est affectée si i ∈ error_rows OU j ∈ error_cols
    
    std::set<int> er_set(error_rows.begin(), error_rows.end());
    std::set<int> ec_set(error_cols.begin(), error_cols.end());
    
    for (int i = 0; i < m; i++) {
        bool row_is_free = (er_set.count(i) > 0);
        
        for (int j = 0; j < n; j++) {
            if (M(i, j) == -1) continue;  // Don't care
            
            bool col_is_free = (ec_set.count(j) > 0);
            
            // Si ni la ligne ni la colonne n'est libre, skip (valeur fixée)
            if (!row_is_free && !col_is_free) continue;
            
            bool target_value = (M(i, j) == 1);
            bool current_value = (computed(i, j) == 1);
            bool is_error = (target_value != current_value);
            
            if (target_value) {
                // M[i,j] = 1 : besoin d'au moins un facteur l tel que A[i,l]=1 ET B[l,j]=1
                std::vector<int> big_or;
                
                for (int l = 0; l < k; l++) {
                    bool a_free = row_is_free;
                    bool b_free = col_is_free;
                    
                    if (!a_free && !b_free) {
                        // Les deux fixés - contribution constante
                        if (A_current(i, l) == 1 && B_current(l, j) == 1) {
                            // Clause toujours satisfaite
                            big_or.clear();
                            big_or.push_back(ultra_solver.new_var());
                            ultra_solver.add_clause({big_or[0]});  // Variable vraie
                            break;
                        }
                    } else if (!a_free) {
                        // A fixé, B libre
                        if (A_current(i, l) == 1) {
                            big_or.push_back(B_vars.get(l, j));
                        }
                    } else if (!b_free) {
                        // A libre, B fixé
                        if (B_current(l, j) == 1) {
                            big_or.push_back(A_vars.get(i, l));
                        }
                    } else {
                        // Les deux libres - variable auxiliaire pour A[i,l] ∧ B[l,j]
                        int var_a = A_vars.get(i, l);
                        int var_b = B_vars.get(l, j);
                        int aux = ultra_solver.new_var();
                        
                        // aux <=> (var_a ∧ var_b)
                        ultra_solver.add_clause({-aux, var_a});
                        ultra_solver.add_clause({-aux, var_b});
                        ultra_solver.add_clause({-var_a, -var_b, aux});
                        hard_count += 3;
                        
                        big_or.push_back(aux);
                    }
                }
                
                if (!big_or.empty()) {
                    if (is_error) {
                        // C'est une erreur actuelle - SOFT clause
                        ultra_solver.add_soft_clause(big_or, 1);
                        soft_count++;
                    } else {
                        // Pas une erreur - HARD clause (préserver)
                        ultra_solver.add_clause(big_or);
                        hard_count++;
                    }
                }
                
            } else {
                // M[i,j] = 0 : pour tout l, NOT(A[i,l] ∧ B[l,j])
                for (int l = 0; l < k; l++) {
                    bool a_free = row_is_free;
                    bool b_free = col_is_free;
                    
                    std::vector<int> clause;
                    
                    if (!a_free && !b_free) {
                        // Les deux fixés - vérifier cohérence
                        if (A_current(i, l) == 1 && B_current(l, j) == 1) {
                            // Erreur fixée - ne peut pas être corrigée ici
                        }
                        continue;
                    } else if (!a_free) {
                        if (A_current(i, l) == 1) {
                            clause = {-B_vars.get(l, j)};
                        }
                    } else if (!b_free) {
                        if (B_current(l, j) == 1) {
                            clause = {-A_vars.get(i, l)};
                        }
                    } else {
                        clause = {-A_vars.get(i, l), -B_vars.get(l, j)};
                    }
                    
                    if (!clause.empty()) {
                        if (is_error) {
                            ultra_solver.add_soft_clause(clause, 1);
                            soft_count++;
                        } else {
                            ultra_solver.add_clause(clause);
                            hard_count++;
                        }
                    }
                }
            }
        }
    }
    
    std::cout << "  [ULTRA-CIBLÉ] Hard: " << hard_count << ", Soft: " << soft_count << std::endl;
    
    // 4. Fixer les variables A et B qui ne sont pas libres mais ont été créées
    std::vector<int> assumptions;
    
    // Fixer A pour les lignes NON dans error_rows
    for (auto& [key, var] : A_vars.vars) {
        int i = key.first;
        int l = key.second;
        if (er_set.count(i) == 0) {
            assumptions.push_back(A_current(i, l) == 1 ? var : -var);
        }
    }
    
    // Fixer B pour les colonnes NON dans error_cols
    for (auto& [key, var] : B_vars.vars) {
        int l = key.first;
        int j = key.second;
        if (ec_set.count(j) == 0) {
            assumptions.push_back(B_current(l, j) == 1 ? var : -var);
        }
    }
    
    // 5. Résoudre
    if (!ultra_solver.solve_with_assumptions(assumptions)) {
        std::cout << "  [ULTRA-CIBLÉ] UNSAT !" << std::endl;
        return -1;
    }
    
    long long cost = ultra_solver.get_cost();
    
    // 6. Extraire la solution
    for (int i : error_rows) {
        for (int l = 0; l < k; l++) {
            auto it = A_vars.vars.find({i, l});
            if (it != A_vars.vars.end()) {
                A_current(i, l) = ultra_solver.getValue(it->second) ? 1 : 0;
            }
        }
    }
    
    for (int l = 0; l < k; l++) {
        for (int j : error_cols) {
            auto it = B_vars.vars.find({l, j});
            if (it != B_vars.vars.end()) {
                B_current(l, j) = ultra_solver.getValue(it->second) ? 1 : 0;
            }
        }
    }
    
    std::cout << "  [ULTRA-CIBLÉ] Résolu avec coût: " << cost << std::endl;
    return static_cast<int>(cost);
}

// ==================== MODE PARTIAL ====================

int BMF::lns_step_partial(
    Matrix& A_current, Matrix& B_current,
    const std::vector<int>& rows_to_free,
    const std::vector<int>& cols_to_free
) {
    // Construire les vecteurs de lignes/colonnes libres
    std::vector<bool> free_rows(m, false);
    std::vector<bool> free_cols(n, false);
    
    for (int r : rows_to_free) if (r >= 0 && r < m) free_rows[r] = true;
    for (int c : cols_to_free) if (c >= 0 && c < n) free_cols[c] = true;
    
    // 1. Calculer le coût fixe (erreurs dans les cellules NON impactées)
    // Une cellule est NON impactée si sa ligne ET sa colonne sont fixées
    int fixed_cost = 0;
    Matrix computed = A_current.multiply(B_current);
    
    for (int i = 0; i < m; i++) {
        if (free_rows[i]) continue;  // Ligne libérée = impactée
        for (int j = 0; j < n; j++) {
            if (free_cols[j]) continue;  // Colonne libérée = impactée
            if (M(i, j) == -1) continue;  // Ignorer les valeurs manquantes (don't care)
            // Cellule (i,j) est complètement fixée
            if (computed(i, j) != M(i, j)) {
                fixed_cost++;
            }
        }
    }
    
    // 2. Créer un nouveau solveur  pour les cellules impactées uniquement
    SATSolver partial_solver;
    VariableManager A_partial(partial_solver), B_partial(partial_solver);
    
    // Créer les variables pour les lignes/colonnes libres
    for (int i = 0; i < m; i++) {
        if (free_rows[i]) {
            for (int l = 0; l < k; l++) {
                A_partial.get(i, l);
            }
        }
    }
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            if (free_cols[j]) {
                B_partial.get(l, j);
            }
        }
    }
    
    int hard_count = 0, soft_count = 0;
    
    // 3. Pour chaque cellule IMPACTÉE, créer les clauses appropriées
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            bool row_free = free_rows[i];
            bool col_free = free_cols[j];
            
            // Cellule non impactée - déjà comptée dans fixed_cost
            if (!row_free && !col_free) continue;
            
            // Ignorer les valeurs manquantes (don't care) - pas de contrainte
            if (M(i, j) == -1) continue;
            
            bool is_error = (computed(i, j) != M(i, j));
            
            if (M(i, j) == 1) {
                // Pour M[i,j] = 1: OR_l (A[i,l] AND B[l,j]) = 1
                std::vector<int> big_or;
                
                for (int l = 0; l < k; l++) {
                    if (!row_free) {
                        // A[i,l] est fixé
                        if (A_current(i, l) == 1) {
                            int var_b = B_partial.get(l, j);
                            big_or.push_back(var_b);
                        }
                    } else if (!col_free) {
                        // B[l,j] est fixé
                        if (B_current(l, j) == 1) {
                            int var_a = A_partial.get(i, l);
                            big_or.push_back(var_a);
                        }
                    } else {
                        // Les deux sont libres - variable auxiliaire
                        int var_a = A_partial.get(i, l);
                        int var_b = B_partial.get(l, j);
                        int aux = partial_solver.new_var();
                        
                        partial_solver.add_clause({-aux, var_a});
                        partial_solver.add_clause({-aux, var_b});
                        partial_solver.add_clause({-var_a, -var_b, aux});
                        hard_count += 3;
                        
                        big_or.push_back(aux);
                    }
                }
                
                if (big_or.empty()) {
                    // Clause impossible à satisfaire
                    if (is_error) {
                        // Erreur inévitable, comptée comme soft
                        int dummy = partial_solver.new_var();
                        partial_solver.add_clause({dummy});
                        partial_solver.add_soft_clause({-dummy}, 1);
                        soft_count++;
                    }
                } else {
                    if (is_error) {
                        partial_solver.add_soft_clause(big_or, 1);
                        soft_count++;
                    } else {
                        partial_solver.add_clause(big_or);
                        hard_count++;
                    }
                }
                
            } else { // M(i, j) == 0
                // Pour M[i,j] = 0: AND_l NOT(A[i,l] AND B[l,j])
                for (int l = 0; l < k; l++) {
                    std::vector<int> clause;
                    
                    if (!row_free) {
                        if (A_current(i, l) == 1) {
                            int var_b = B_partial.get(l, j);
                            clause = {-var_b};
                        }
                        // Si A[i,l] = 0, clause toujours satisfaite
                    } else if (!col_free) {
                        if (B_current(l, j) == 1) {
                            int var_a = A_partial.get(i, l);
                            clause = {-var_a};
                        }
                    } else {
                        int var_a = A_partial.get(i, l);
                        int var_b = B_partial.get(l, j);
                        clause = {-var_a, -var_b};
                    }
                    
                    if (!clause.empty()) {
                        if (is_error) {
                            partial_solver.add_soft_clause(clause, 1);
                            soft_count++;
                        } else {
                            partial_solver.add_clause(clause);
                            hard_count++;
                        }
                    }
                }
            }
        }
    }
    
    std::cout << "  [PARTIAL] Fixed cost: " << fixed_cost 
              << ", Hard: " << hard_count << ", Soft: " << soft_count << std::endl;
    
    // 4. Ajouter les assumptions pour fixer les variables des lignes/colonnes non-libres
    // (mais qui ont des variables car ils partagent une colonne/ligne avec une libérée)
    std::vector<int> assumptions;
    
    // Fixer A pour les lignes non-libres (mais qui ont des colonnes libres)
    for (int i = 0; i < m; i++) {
        if (free_rows[i]) continue;
        for (int l = 0; l < k; l++) {
            auto it = A_partial.vars.find({i, l});
            if (it != A_partial.vars.end()) {
                assumptions.push_back(A_current(i, l) == 1 ? it->second : -it->second);
            }
        }
    }
    
    // Fixer B pour les colonnes non-libres
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            if (free_cols[j]) continue;
            auto it = B_partial.vars.find({l, j});
            if (it != B_partial.vars.end()) {
                assumptions.push_back(B_current(l, j) == 1 ? it->second : -it->second);
            }
        }
    }
    
    // 5. Résoudre
    if (!partial_solver.solve_with_assumptions(assumptions)) {
        return -1;
    }
    
    long long partial_cost = partial_solver.get_cost();
    
    // 6. Extraire la solution pour les variables libres
    for (int i = 0; i < m; i++) {
        if (!free_rows[i]) continue;
        for (int l = 0; l < k; l++) {
            auto it = A_partial.vars.find({i, l});
            if (it != A_partial.vars.end()) {
                A_current(i, l) = partial_solver.getValue(it->second) ? 1 : 0;
            }
        }
    }
    
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            if (!free_cols[j]) continue;
            auto it = B_partial.vars.find({l, j});
            if (it != B_partial.vars.end()) {
                B_current(l, j) = partial_solver.getValue(it->second) ? 1 : 0;
            }
        }
    }
    
    return fixed_cost + static_cast<int>(partial_cost);
}


LNSResult BMF::solve_lns_v3(
    int max_iterations,
    double neighborhood_size,
    unsigned seed,
    int stagnation_threshold,
    bool verbose,
    double timeout_ms,
    std::atomic<bool>* stop_flag
) {
    LNSResult result;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Lambda pour vérifier les conditions d'arrêt
    auto should_stop = [&]() -> bool {
        // 1. Vérifier le flag d'interruption externe (Ctrl+C)
        if (stop_flag && stop_flag->load()) {
            result.stop_reason = "interrupted";
            return true;
        }
        // 2. Vérifier le timeout
        if (timeout_ms > 0) {
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double, std::milli>(now - start_time).count();
            if (elapsed >= timeout_ms) {
                result.stop_reason = "timeout";
                return true;
            }
        }
        return false;
    };
    
    std::mt19937 rng(seed);
    reset();
    
    // === 1. INITIALISATION ===
    Matrix A_current(m, k, 0);
    Matrix B_current(k, n, 0);
    double greedy_time = 0;
    
    // Vérifier si une solution initiale a été fournie via set_initial_solution
    if (has_initial_solution && initial_A.rows > 0 && initial_B.rows > 0) {
        // Utiliser la solution fournie
        A_current = initial_A;
        B_current = initial_B;
        result.init_method = "from_solution+lns_v3";
        if (verbose) {
            std::cout << "Using provided initial solution for LNS-v3." << std::endl;
        }
    } else {
        // Initialisation gloutonne via BMFLocalSearch
        auto greedy_start = std::chrono::high_resolution_clock::now();
        BMFLocalSearch ls(k, M, seed);
        ls.initialize_greedy();
        auto greedy_end = std::chrono::high_resolution_clock::now();
        greedy_time = std::chrono::duration<double, std::milli>(greedy_end - greedy_start).count();
        
        A_current = ls.A;
        B_current = ls.B;
        result.init_method = "greedy+lns_v3";
    }
    
    // === 2. ÉVALUATION INITIALE ===
    int best_cost = count_alignment_errors(A_current, B_current);
    Matrix A_best = A_current;
    Matrix B_best = B_current;
    
    result.initial_errors = best_cost;
    result.error_history.push_back(best_cost);
    
    // Calculer la taille du voisinage en avance pour l'affichage
    //int rows_to_free = std::max(2, static_cast<int>(m * neighborhood_size));
    //int cols_to_free = std::max(2, static_cast<int>(n * neighborhood_size));
    int rows_to_free = 15;
    int cols_to_free = 15;
    
    if (verbose) {
        std::cout << "\n=== LNS-v3: Alternance A/B + Anti-stagnation ===" << std::endl;
        std::cout << "Initial errors: " << best_cost << std::endl;
        std::cout << "Greedy initialization time: " << std::fixed << std::setprecision(2) << greedy_time << " ms" << std::endl;
        std::cout << "Neighborhood size: " << rows_to_free << "x" << cols_to_free 
                  << " (" << (neighborhood_size * 100) << "% of " << m << "x" << n << ")" << std::endl;
        std::cout << "Stagnation threshold: " << stagnation_threshold << std::endl;
    }
    
    if (best_cost == 0) {
        auto end_time = std::chrono::high_resolution_clock::now();
        result.total_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        result.success = true;
        result.final_errors = 0;
        result.A_solution = A_best;
        result.B_solution = B_best;
        if (verbose) {
            std::cout << "\n*** FACTORISATION EXACTE trouvée par initialisation gloutonne! ***" << std::endl;
            std::cout << "Temps greedy: " << std::fixed << std::setprecision(2) << greedy_time << " ms" << std::endl;
            std::cout << "Aucune itération LNS nécessaire." << std::endl;
        }
        return result;
    }
    
    // === 3. BOUCLE PRINCIPALE LNS-v3 ===
    int stagnation = 0;
    int total_restarts = 0;
    
    // === ESCALADE PROGRESSIVE ===
    // Tailles de voisinage: 8, 12, 15, 20, 25, 40, 50, 70
    const std::vector<int> NEIGHBORHOOD_SIZES = {20, 25, 40, 50, 70, 100};
    int escalation_level = 0;
    int base_neighborhood_idx = 0;  // Index dans NEIGHBORHOOD_SIZES
    
    // Trouver l'index de départ basé sur la taille initiale
    for (size_t i = 0; i < NEIGHBORHOOD_SIZES.size(); i++) {
        if (NEIGHBORHOOD_SIZES[i] >= std::min(rows_to_free, cols_to_free)) {
            base_neighborhood_idx = static_cast<int>(i);
            break;
        }
    }
    
    // Variables pour la destruction de facteurs
    std::set<int> recently_destroyed;
    int iterations_level5 = 0;  // Compteur pour détecter boucle au niveau 5+
    int total_full_restarts = 0;  // Compteur de RESTART TOTAL
    const int MAX_FULL_RESTARTS = 5;  // Maximum de restarts totaux (augmenté encore)
    
    // === NOUVELLES VARIABLES POUR EXPLORATION LONGUE ===
    int iterations_since_last_improvement = 0;  // Compteur global
    int exploration_mode_iterations = 0;  // Iterations restantes en mode exploration
    const int EXPLORATION_MODE_DURATION = 100;  // Durée du mode exploration après restart
    int last_improvement_cost = best_cost;  // Dernier coût où on a eu une amélioration
    
    // Seuil de stagnation adaptatif (plus élevé au début)
    int adaptive_stagnation_threshold = stagnation_threshold * 2;  // Commence à 30
    const int JOINT_PERIOD = 10;  // Phase JOINT tous les 10 itérations
    
    for (int iter = 1; iter <= max_iterations; iter++) {
        // === VÉRIFICATION DES CONDITIONS D'ARRÊT ===
        if (should_stop()) {
            if (verbose) {
                std::cout << "\nArrêt demandé (" << result.stop_reason << ") à l'itération " << iter << std::endl;
            }
            break;
        }
        
        auto iter_start = std::chrono::high_resolution_clock::now();
        
        // === 3.1 DÉTERMINER LA PHASE ===
        OptimizationPhase phase;
        if (iter % JOINT_PERIOD == 0) {
            phase = OptimizationPhase::JOINT;
        } else if (iter % 2 == 0) {
            phase = OptimizationPhase::FIX_A;  // Optimiser B
        } else {
            phase = OptimizationPhase::FIX_B;  // Optimiser A
        }
        
        // === 3.2 SÉLECTION CAUSALE PAR FACTEUR ===
        // Calculer les erreurs par ligne, colonne ET par facteur
        auto [row_errors, col_errors] = compute_errors(A_current, B_current);
        
        // NOUVEAU: Calculer les erreurs causées par chaque facteur
        std::vector<int> factor_errors_under(k, 0);  // Facteurs absents où on a besoin de couverture
        std::vector<int> factor_errors_over(k, 0);   // Facteurs qui causent des overcouvertures
        
        Matrix computed = A_current.multiply(B_current);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (M(i, j) == -1) continue;
                
                if (M(i, j) == 1 && computed(i, j) == 0) {
                    // Undercover: quels facteurs auraient pu couvrir?
                    for (int l = 0; l < k; l++) {
                        // Un facteur peut couvrir si A[i,l] OU B[l,j] est à 1
                        if (A_current(i, l) == 1 || B_current(l, j) == 1) {
                            factor_errors_under[l]++;
                        }
                    }
                } else if (M(i, j) == 0 && computed(i, j) == 1) {
                    // Overcover: quels facteurs causent cette erreur?
                    for (int l = 0; l < k; l++) {
                        if (A_current(i, l) == 1 && B_current(l, j) == 1) {
                            factor_errors_over[l]++;
                        }
                    }
                }
            }
        }
        
        // Identifier les facteurs les plus problématiques
        std::vector<std::pair<int, int>> factor_scores(k);
        for (int l = 0; l < k; l++) {
            // Score = erreurs totales causées ou non couvertes par ce facteur
            factor_scores[l] = {factor_errors_under[l] + factor_errors_over[l], l};
        }
        std::sort(factor_scores.begin(), factor_scores.end(), std::greater<>());
        
        // Les top-2 facteurs problématiques
        std::vector<int> worst_factors;
        for (int f = 0; f < std::min(2, k); f++) {
            if (factor_scores[f].first > 0) {
                worst_factors.push_back(factor_scores[f].second);
            }
        }
        
        std::vector<int> selected_rows, selected_cols;
        
        // === MODE SÉLECTION CIBLÉE (niveau 5+) ===
        // Si on est au niveau 5+, sélectionner uniquement les lignes/colonnes d'erreurs
        bool use_targeted_selection = (escalation_level >= 5 && escalation_level <= 8 && best_cost <= 20);
        
        if (use_targeted_selection) {
            // Résolution ciblée : libérer uniquement les lignes/colonnes impliquées dans les erreurs
            Matrix computed = A_best.multiply(B_best);
            std::set<int> error_rows_set, error_cols_set;
            
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    if (M(i, j) != -1 && M(i, j) != computed(i, j)) {
                        error_rows_set.insert(i);
                        error_cols_set.insert(j);
                    }
                }
            }
            
            // Étendre pour avoir plus de contexte (plus on stagne, plus on étend)
            int extend = 3 + iterations_level5 * 5;  // 8, 13, 18, 23, 28, 33
            
            std::set<int> rows_set, cols_set;
            for (int r : error_rows_set) {
                for (int dr = -extend; dr <= extend; dr++) {
                    int nr = r + dr;
                    if (nr >= 0 && nr < m) rows_set.insert(nr);
                }
            }
            for (int c : error_cols_set) {
                for (int dc = -extend; dc <= extend; dc++) {
                    int nc = c + dc;
                    if (nc >= 0 && nc < n) cols_set.insert(nc);
                }
            }
            
            for (int r : rows_set) selected_rows.push_back(r);
            for (int c : cols_set) selected_cols.push_back(c);
            
        } else {
            // HYBRID: 70% Top-K + 30% random (adaptatif avec stagnation)
            std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
            double random_ratio = 0.3 + 0.2 * (static_cast<double>(stagnation) / stagnation_threshold);
            random_ratio = std::min(random_ratio, 0.6);  // Max 60% random
        
        // Sélection selon la phase
        switch (phase) {
            case OptimizationPhase::FIX_A: {
                // Fixer A, libérer colonnes de B
                // NOUVEAU: Score hybride = erreurs + impact sur pires facteurs
                std::vector<std::pair<int, int>> col_scores;
                for (int j = 0; j < n; j++) {
                    int score = col_errors[j];
                    // Bonus si cette colonne touche aux pires facteurs
                    for (int wf : worst_factors) {
                        if (B_current(wf, j) == 1 && factor_errors_over[wf] > 0) {
                            score += factor_errors_over[wf] / 2;  // Ce facteur cause des overcouvertures
                        }
                        if (B_current(wf, j) == 0 && factor_errors_under[wf] > 0) {
                            score += factor_errors_under[wf] / 2;  // Ce facteur pourrait aider
                        }
                    }
                    col_scores.push_back({score, j});
                }
                std::sort(col_scores.begin(), col_scores.end(), std::greater<>());
                
                for (int j = 0; j < std::min(cols_to_free, n); j++) {
                    if (prob_dist(rng) < random_ratio) {
                        // Random
                        int rand_col = rng() % n;
                        selected_cols.push_back(rand_col);
                    } else {
                        // Top-K (avec score hybride)
                        selected_cols.push_back(col_scores[j].second);
                    }
                }
                break;
            }
            
            case OptimizationPhase::FIX_B: {
                // Fixer B, libérer lignes de A
                // NOUVEAU: Score hybride = erreurs + impact sur pires facteurs
                std::vector<std::pair<int, int>> row_scores;
                for (int i = 0; i < m; i++) {
                    int score = row_errors[i];
                    // Bonus si cette ligne touche aux pires facteurs
                    for (int wf : worst_factors) {
                        if (A_current(i, wf) == 1 && factor_errors_over[wf] > 0) {
                            score += factor_errors_over[wf] / 2;
                        }
                        if (A_current(i, wf) == 0 && factor_errors_under[wf] > 0) {
                            score += factor_errors_under[wf] / 2;
                        }
                    }
                    row_scores.push_back({score, i});
                }
                std::sort(row_scores.begin(), row_scores.end(), std::greater<>());
                
                for (int i = 0; i < std::min(rows_to_free, m); i++) {
                    if (prob_dist(rng) < random_ratio) {
                        int rand_row = rng() % m;
                        selected_rows.push_back(rand_row);
                    } else {
                        selected_rows.push_back(row_scores[i].second);
                    }
                }
                // selected_cols reste vide (B est fixé entièrement)
                break;
            }
            
            case OptimizationPhase::JOINT: {
                // Libérer un sous-ensemble des deux (voisinage réduit)
                int joint_rows = std::max(1, rows_to_free / 2);
                int joint_cols = std::max(1, cols_to_free / 2);
                
                std::vector<std::pair<int, int>> row_scores, col_scores;
                for (int i = 0; i < m; i++) row_scores.push_back({row_errors[i], i});
                for (int j = 0; j < n; j++) col_scores.push_back({col_errors[j], j});
                std::sort(row_scores.begin(), row_scores.end(), std::greater<>());
                std::sort(col_scores.begin(), col_scores.end(), std::greater<>());
                
                for (int i = 0; i < std::min(joint_rows, m); i++) {
                    selected_rows.push_back(row_scores[i].second);
                }
                for (int j = 0; j < std::min(joint_cols, n); j++) {
                    selected_cols.push_back(col_scores[j].second);
                }
                break;
            }
        }
        }  // Fin du else (mode normal vs ciblé)
        
        // === 3.2.5 PRÉ-RÉPARATION DES LIGNES/COLONNES MORTES ===
        // Avant FIX_A: réparer les lignes de A entièrement à zéro
        if (phase == OptimizationPhase::FIX_A) {
            for (int i = 0; i < m; i++) {
                bool has_active = false;
                for (int l = 0; l < k; l++) {
                    if (A_current(i, l) == 1) { has_active = true; break; }
                }
                
                if (!has_active) {
                    // Ligne i morte - choisir le meilleur facteur à activer
                    int best_l = -1, best_gain = INT_MIN;
                    for (int l = 0; l < k; l++) {
                        int gain = 0;
                        for (int j = 0; j < n; j++) {
                            if (M(i, j) == -1) continue;
                            if (M(i, j) == 1 && B_current(l, j) == 1) gain++;
                            if (M(i, j) == 0 && B_current(l, j) == 1) gain--;
                        }
                        if (gain > best_gain) {
                            best_gain = gain;
                            best_l = l;
                        }
                    }
                    if (best_l >= 0 && best_gain > 0) {
                        A_current(i, best_l) = 1;  // Réparer la ligne morte
                    }
                }
            }
        }
        // Avant FIX_B: réparer les colonnes de B entièrement à zéro
        else if (phase == OptimizationPhase::FIX_B) {
            for (int j = 0; j < n; j++) {
                bool has_active = false;
                for (int l = 0; l < k; l++) {
                    if (B_current(l, j) == 1) { has_active = true; break; }
                }
                
                if (!has_active) {
                    // Colonne j morte - choisir le meilleur facteur à activer
                    int best_l = -1, best_gain = INT_MIN;
                    for (int l = 0; l < k; l++) {
                        int gain = 0;
                        for (int i = 0; i < m; i++) {
                            if (M(i, j) == -1) continue;
                            if (M(i, j) == 1 && A_current(i, l) == 1) gain++;
                            if (M(i, j) == 0 && A_current(i, l) == 1) gain--;
                        }
                        if (gain > best_gain) {
                            best_gain = gain;
                            best_l = l;
                        }
                    }
                    if (best_l >= 0 && best_gain > 0) {
                        B_current(best_l, j) = 1;  // Réparer la colonne morte
                    }
                }
            }
        }
        
        // === 3.3 RÉSOLUTION MaxSAT SELON LA PHASE ===
        int new_cost;
        
        if (phase == OptimizationPhase::FIX_A) {
            // Optimiser B avec A fixé
            // Créer un solveur partiel pour B seulement
            new_cost = lns_step_fix_A(A_current, B_current, selected_cols);
        } else if (phase == OptimizationPhase::FIX_B) {
            // Optimiser A avec B fixé
            new_cost = lns_step_fix_B(A_current, B_current, selected_rows);
        } else {
            // JOINT: utiliser le mode PARTIAL existant
            new_cost = lns_step_partial(A_current, B_current, selected_rows, selected_cols);
        }
        
        auto iter_end = std::chrono::high_resolution_clock::now();
        double iter_time = std::chrono::duration<double, std::milli>(iter_end - iter_start).count();
        result.iteration_times.push_back(iter_time);
        
        // === 3.4 MISE À JOUR ===
        if (new_cost >= 0 && new_cost < best_cost) {
            if (verbose) {
                std::string phase_name = (phase == OptimizationPhase::FIX_A) ? "FIX_A" :
                                         (phase == OptimizationPhase::FIX_B) ? "FIX_B" : "JOINT";
                std::cout << "Iter " << iter << " [" << phase_name << "]: " 
                          << best_cost << " -> " << new_cost << std::endl;
            }
            best_cost = new_cost;
            A_best = A_current;
            B_best = B_current;
            stagnation = 0;
            iterations_since_last_improvement = 0;
            
            // Activer le mode exploration si amélioration significative
            if (best_cost < last_improvement_cost) {
                last_improvement_cost = best_cost;
                exploration_mode_iterations = EXPLORATION_MODE_DURATION;
                // Augmenter le seuil de stagnation pour explorer plus
                adaptive_stagnation_threshold = stagnation_threshold * 2;
            }
            
            // Reset de l'escalade après amélioration
            escalation_level = 0;
            recently_destroyed.clear();
            iterations_level5 = 0;  // Reset du compteur de boucle
            
            if (best_cost == 0) {
                result.stop_reason = "zero_errors";
                if (verbose) std::cout << "Solution parfaite trouvée!" << std::endl;
                break;
            }
        } else {
            stagnation++;
            iterations_since_last_improvement++;
            
            // Réduire le mode exploration
            if (exploration_mode_iterations > 0) {
                exploration_mode_iterations--;
            }
            
            // Adapter le seuil de stagnation progressivement
            if (iterations_since_last_improvement > 200 && adaptive_stagnation_threshold > stagnation_threshold) {
                adaptive_stagnation_threshold = std::max(stagnation_threshold, adaptive_stagnation_threshold - 1);
            }
            
            // === 3.5 ANTI-STAGNATION AVEC ESCALADE PROGRESSIVE ===
            // Utiliser le seuil adaptatif (ou le seuil normal si en mode exploration)
            int effective_threshold = (exploration_mode_iterations > 0) ? 
                                      stagnation_threshold * 3 : adaptive_stagnation_threshold;
            
            if (stagnation >= effective_threshold) {
                total_restarts++;
                escalation_level++;
                
                // === STRATÉGIE D'ESCALADE ===
                // Niveau 1-2: Voisinage élargi progressif
                // Niveau 3-4: DESTRUCTION-RECONSTRUCTION de facteur
                // Niveau 5+: Greedy résidu + voisinage maximal
                
                if (escalation_level <= 2) {
                    // === NIVEAU 1-2: ESCALADE PROGRESSIVE DU VOISINAGE ===
                    int new_idx = std::min(base_neighborhood_idx + escalation_level, 
                                          static_cast<int>(NEIGHBORHOOD_SIZES.size()) - 1);
                    int new_size = NEIGHBORHOOD_SIZES[new_idx];
                    
                    rows_to_free = std::min(new_size, m);
                    cols_to_free = std::min(new_size, n);
                    
                    if (verbose) {
                        std::cout << "Iter " << iter << ": ESCALADE niveau " << escalation_level
                                  << " - voisinage " << rows_to_free << "x" << cols_to_free << std::endl;
                    }
                    
                    // Perturbation légère centrée sur les erreurs
                    double perturb_rate = 0.10 + 0.05 * escalation_level;
                    Matrix current_computed = A_best.multiply(B_best);
                    
                    std::uniform_real_distribution<double> flip_dist(0.0, 1.0);
                    
                    // Pour chaque ligne avec erreurs, trouver le pire facteur
                    for (int i = 0; i < m; i++) {
                        if (row_errors[i] > 0) {
                            // Calculer l'impact de chaque facteur sur cette ligne
                            std::vector<int> factor_impact(k, 0);
                            for (int l = 0; l < k; l++) {
                                for (int j = 0; j < n; j++) {
                                    if (M(i, j) == -1) continue;
                                    
                                    bool currently_covered = (current_computed(i, j) == 1);
                                    
                                    // Que se passe-t-il si on flip A[i,l] ?
                                    bool would_cover = false;
                                    int new_A_il = 1 - A_best(i, l);
                                    for (int ll = 0; ll < k; ll++) {
                                        int A_il_val = (ll == l) ? new_A_il : A_best(i, ll);
                                        if (A_il_val == 1 && B_best(ll, j) == 1) {
                                            would_cover = true;
                                            break;
                                        }
                                    }
                                    
                                    // Impact sur l'erreur
                                    if (M(i, j) == 1) {  // Devrait être couvert
                                        if (!currently_covered && would_cover) factor_impact[l]++;
                                        if (currently_covered && !would_cover) factor_impact[l]--;
                                    } else {  // Ne devrait pas être couvert
                                        if (currently_covered && !would_cover) factor_impact[l]++;
                                        if (!currently_covered && would_cover) factor_impact[l]--;
                                    }
                                }
                            }
                            
                            // Trouver le meilleur facteur à flipper
                            int best_l = 0;
                            int best_impact = factor_impact[0];
                            for (int l = 1; l < k; l++) {
                                if (factor_impact[l] > best_impact) {
                                    best_impact = factor_impact[l];
                                    best_l = l;
                                }
                            }
                            
                            // Flipper le meilleur facteur avec haute probabilité
                            for (int l = 0; l < k; l++) {
                                if (l == best_l && best_impact > 0) {
                                    // Flip ciblé: haute probabilité
                                    if (flip_dist(rng) < 0.8) {
                                        A_current(i, l) = 1 - A_best(i, l);
                                    } else {
                                        A_current(i, l) = A_best(i, l);
                                    }
                                } else if (flip_dist(rng) < perturb_rate) {
                                    // Flip aléatoire pour exploration
                                    A_current(i, l) = 1 - A_best(i, l);
                                } else {
                                    A_current(i, l) = A_best(i, l);
                                }
                            }
                        } else {
                            for (int l = 0; l < k; l++) {
                                A_current(i, l) = A_best(i, l);
                            }
                        }
                    }
                    
                    // Pour chaque colonne avec erreurs, trouver le pire facteur
                    for (int j = 0; j < n; j++) {
                        if (col_errors[j] > 0) {
                            std::vector<int> factor_impact(k, 0);
                            for (int l = 0; l < k; l++) {
                                for (int i = 0; i < m; i++) {
                                    if (M(i, j) == -1) continue;
                                    
                                    bool currently_covered = (current_computed(i, j) == 1);
                                    
                                    // Que se passe-t-il si on flip B[l,j] ?
                                    bool would_cover = false;
                                    int new_B_lj = 1 - B_best(l, j);
                                    for (int ll = 0; ll < k; ll++) {
                                        int B_lj_val = (ll == l) ? new_B_lj : B_best(ll, j);
                                        if (A_best(i, ll) == 1 && B_lj_val == 1) {
                                            would_cover = true;
                                            break;
                                        }
                                    }
                                    
                                    if (M(i, j) == 1) {
                                        if (!currently_covered && would_cover) factor_impact[l]++;
                                        if (currently_covered && !would_cover) factor_impact[l]--;
                                    } else {
                                        if (currently_covered && !would_cover) factor_impact[l]++;
                                        if (!currently_covered && would_cover) factor_impact[l]--;
                                    }
                                }
                            }
                            
                            int best_l = 0;
                            int best_impact = factor_impact[0];
                            for (int l = 1; l < k; l++) {
                                if (factor_impact[l] > best_impact) {
                                    best_impact = factor_impact[l];
                                    best_l = l;
                                }
                            }
                            
                            for (int l = 0; l < k; l++) {
                                if (l == best_l && best_impact > 0) {
                                    if (flip_dist(rng) < 0.8) {
                                        B_current(l, j) = 1 - B_best(l, j);
                                    } else {
                                        B_current(l, j) = B_best(l, j);
                                    }
                                } else if (flip_dist(rng) < perturb_rate) {
                                    B_current(l, j) = 1 - B_best(l, j);
                                } else {
                                    B_current(l, j) = B_best(l, j);
                                }
                            }
                        } else {
                            for (int l = 0; l < k; l++) {
                                B_current(l, j) = B_best(l, j);
                            }
                        }
                    }
                    
                } else if (escalation_level <= 4) {
                    // === NIVEAU 3-4: DESTRUCTION MULTI-FACTEURS SIMULTANÉE ===
                    // Détruire TOUS les facteurs liés aux erreurs EN MÊME TEMPS
                    // C'est la clé pour sortir des minimums locaux
                    
                    // Calculer la pertinence de chaque facteur par rapport aux erreurs
                    Matrix current_computed = A_best.multiply(B_best);
                    std::vector<int> factor_relevance(k, 0);
                    int total_errors = 0;
                    
                    for (int i = 0; i < m; i++) {
                        for (int j = 0; j < n; j++) {
                            if (M(i, j) == -1) continue;
                            bool is_error = (M(i, j) != current_computed(i, j));
                            if (is_error) {
                                total_errors++;
                                for (int l = 0; l < k; l++) {
                                    if (A_best(i, l) == 1 || B_best(l, j) == 1) {
                                        factor_relevance[l]++;
                                    }
                                }
                            }
                        }
                    }
                    
                    // Collecter TOUS les facteurs avec relevance > 0
                    std::vector<int> factors_to_destroy;
                    for (int l = 0; l < k; l++) {
                        if (factor_relevance[l] > 0 && recently_destroyed.find(l) == recently_destroyed.end()) {
                            factors_to_destroy.push_back(l);
                        }
                    }
                    
                    // Si tous déjà détruits, reset et prendre les top 5
                    if (factors_to_destroy.empty()) {
                        recently_destroyed.clear();
                        std::vector<std::pair<int, int>> sorted_factors;
                        for (int l = 0; l < k; l++) {
                            if (factor_relevance[l] > 0) {
                                sorted_factors.push_back({factor_relevance[l], l});
                            }
                        }
                        std::sort(sorted_factors.rbegin(), sorted_factors.rend());
                        for (size_t i = 0; i < std::min(sorted_factors.size(), size_t(5)); i++) {
                            factors_to_destroy.push_back(sorted_factors[i].second);
                        }
                    }
                    
                    // Limiter à max 2 facteurs pour que MaxSAT reste rapide
                    // (5 facteurs × 150×150 = problèmes avec 1.6M+ clauses hard)
                    if (factors_to_destroy.size() > 2) {
                        // Garder seulement les 2 plus pertinents
                        std::vector<std::pair<int, int>> sorted;
                        for (int l : factors_to_destroy) {
                            sorted.push_back({factor_relevance[l], l});
                        }
                        std::sort(sorted.rbegin(), sorted.rend());
                        factors_to_destroy.clear();
                        for (int i = 0; i < 2 && i < (int)sorted.size(); i++) {
                            factors_to_destroy.push_back(sorted[i].second);
                        }
                    }
                    
                    if (verbose) {
                        std::cout << "Iter " << iter << ": DESTRUCTION LÉGÈRE [";
                        for (size_t i = 0; i < factors_to_destroy.size(); i++) {
                            if (i > 0) std::cout << ",";
                            std::cout << factors_to_destroy[i];
                        }
                        std::cout << "] (" << factors_to_destroy.size() << " facteurs, " 
                                  << total_errors << " erreurs)" << std::endl;
                    }
                    
                    // Partir de la meilleure solution
                    A_current = A_best;
                    B_current = B_best;
                    
                    // DESTRUCTION PARTIELLE: perturber au lieu de mettre à 0
                    // Mettre à 0 seulement 30% des cellules de chaque facteur (au hasard)
                    for (int l : factors_to_destroy) {
                        for (int i = 0; i < m; i++) {
                            if (rng() % 100 < 30) A_current(i, l) = 0;
                        }
                        for (int j = 0; j < n; j++) {
                            if (rng() % 100 < 30) B_current(l, j) = 0;
                        }
                        recently_destroyed.insert(l);
                    }
                    
                    // Voisinage modéré pour reconstruction (pas 150)
                    rows_to_free = std::min(60, m);
                    cols_to_free = std::min(60, n);
                    
                } else if (escalation_level <= 8) {
                    // === NIVEAU 5-8: RÉSOLUTION ULTRA-CIBLÉE ===
                    // Utiliser la méthode ultra-ciblée qui libère UNIQUEMENT
                    // les lignes/colonnes impliquées dans les erreurs
                    iterations_level5++;
                    
                    // Détection de boucle rapide: 3 itérations max au niveau 5+
                    if (iterations_level5 > 3) {
                        if (verbose) {
                            std::cout << "Iter " << iter << ": Stagnation niveau 5+ (" 
                                      << iterations_level5 << " iter) -> RESTART" << std::endl;
                        }
                        escalation_level = 9;  // Forcer RESTART TOTAL
                        iterations_level5 = 0;
                    } else if (best_cost <= 15) {
                        // Peu d'erreurs: utiliser la résolution ultra-ciblée
                        if (verbose) {
                            std::cout << "Iter " << iter << ": RÉSOLUTION ULTRA-CIBLÉE "
                                      << "(niveau 5+, " << iterations_level5 << "/3, " 
                                      << best_cost << " erreurs)" << std::endl;
                        }
                        
                        // Appeler directement la résolution ultra-ciblée
                        A_current = A_best;
                        B_current = B_best;
                        int ultra_cost = lns_step_ultra_targeted(A_current, B_current);
                        
                        if (ultra_cost >= 0 && ultra_cost < best_cost) {
                            if (verbose) {
                                std::cout << "  Amélioration ultra-ciblée: " << best_cost 
                                          << " -> " << ultra_cost << std::endl;
                            }
                            best_cost = ultra_cost;
                            A_best = A_current;
                            B_best = B_current;
                            stagnation = 0;
                            escalation_level = 0;
                            recently_destroyed.clear();
                            iterations_level5 = 0;
                            
                            if (best_cost == 0) {
                                result.stop_reason = "zero_errors";
                                if (verbose) std::cout << "Solution parfaite trouvée!" << std::endl;
                                break;
                            }
                        }
                        // Skip le reste de l'itération (on a déjà fait la résolution)
                        result.error_history.push_back(best_cost);
                        continue;
                    } else {
                        // Plus de 15 erreurs: voisinage large classique
                        rows_to_free = std::min(100, m);
                        cols_to_free = std::min(100, n);
                        if (verbose) {
                            std::cout << "Iter " << iter << ": Voisinage large " 
                                      << rows_to_free << "x" << cols_to_free
                                      << " (" << best_cost << " erreurs)" << std::endl;
                        }
                    }
                } else {
                    // === NIVEAU 9+: RESTART TOTAL ===
                    total_full_restarts++;
                    
                    // Au lieu d'arrêter, on continue avec des stratégies plus agressives
                    // après MAX_FULL_RESTARTS, on essaie des perturbations aléatoires plus fortes
                    bool beyond_normal_restarts = (total_full_restarts > MAX_FULL_RESTARTS);
                    
                    if (verbose) {
                        if (beyond_normal_restarts) {
                            std::cout << "Iter " << iter << ": RESTART AGRESSIF " << total_full_restarts 
                                      << " (au-delà de " << MAX_FULL_RESTARTS << ", erreurs: " << best_cost << ")" << std::endl;
                        } else {
                            std::cout << "Iter " << iter << ": RESTART TOTAL " << total_full_restarts 
                                      << "/" << MAX_FULL_RESTARTS << " (erreurs: " << best_cost << ")" << std::endl;
                        }
                    }
                    
                    // Nouvelle seed basée sur l'itération et le numéro de restart
                    std::mt19937 restart_rng(seed + iter * 1000 + total_full_restarts * 7919);
                    std::uniform_int_distribution<int> dist(0, 1);
                    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
                    
                    // Identifier les facteurs impliqués dans les erreurs (utilisé par toutes les stratégies)
                    Matrix computed = A_best.multiply(B_best);
                    std::set<int> error_factors;
                    std::set<int> error_rows_set, error_cols_set;
                    std::vector<std::pair<int, int>> errors;  // Liste des cellules en erreur
                    for (int i = 0; i < m; i++) {
                        for (int j = 0; j < n; j++) {
                            if (M(i, j) != -1 && M(i, j) != computed(i, j)) {
                                error_rows_set.insert(i);
                                error_cols_set.insert(j);
                                errors.push_back({i, j});
                                for (int l = 0; l < k; l++) {
                                    if (A_best(i, l) == 1 || B_best(l, j) == 1) {
                                        error_factors.insert(l);
                                    }
                                }
                            }
                        }
                    }
                    
                    // Stratégie différente selon le numéro de restart
                    // NOUVELLE ROTATION OPTIMISÉE: S3 et S6 sont les plus efficaces (faibles erreurs)
                    // Rotation: S3->S6->S3->S6->S8->S3->S6->S1 (S8 = nouvelle stratégie intelligente)
                    // Les stratégies S0, S2, S4, S5, S7 créent trop d'erreurs
                    static const int strategy_order[] = {3, 6, 3, 6, 8, 3, 6, 1};
                    int strategy = strategy_order[total_full_restarts % 8];
                    
                    // Limiter le nombre de facteurs perturbés (max 3 pour garder le problème simple)
                    std::vector<int> limited_error_factors;
                    for (int l : error_factors) {
                        limited_error_factors.push_back(l);
                        if (limited_error_factors.size() >= 3) break;  // Max 3 facteurs
                    }
                    
                    if (strategy == 8) {
                        // === NOUVELLE STRATÉGIE S8: Analyse intelligente des patterns d'erreur ===
                        // Idée: Identifier les cellules problématiques et modifier UNIQUEMENT
                        // les facteurs spécifiques qui causent chaque erreur
                        A_current = A_best;
                        B_current = B_best;
                        
                        // Diagnostic: afficher les 5 erreurs restantes
                        if (verbose && errors.size() <= 10) {
                            std::cout << "  [DIAGNOSTIC] " << errors.size() << " erreurs à corriger:" << std::endl;
                            for (auto& [ei, ej] : errors) {
                                bool expected = (M(ei, ej) == 1);
                                bool got = (computed(ei, ej) == 1);
                                std::cout << "    - Cellule (" << ei << "," << ej << "): attendu=" 
                                          << expected << ", obtenu=" << got;
                                // Compter les facteurs actifs pour cette cellule
                                int active_factors = 0;
                                std::vector<int> active_factor_list;
                                for (int l = 0; l < k; l++) {
                                    if (A_best(ei, l) == 1 && B_best(l, ej) == 1) {
                                        active_factors++;
                                        if (active_factor_list.size() < 3) active_factor_list.push_back(l);
                                    }
                                }
                                std::cout << " (" << active_factors << " facteurs actifs)" << std::endl;
                            }
                        }
                        
                        // Pour chaque erreur, appliquer une correction ciblée
                        for (auto& [ei, ej] : errors) {
                            bool expected = (M(ei, ej) == 1);
                            
                            if (expected) {
                                // Erreur: devrait être 1, mais c'est 0
                                // Solution: activer UN facteur pour cette cellule
                                // Choisir le facteur avec le plus de 1 dans la ligne ei
                                int best_factor = -1;
                                int best_score = -1;
                                for (int l : limited_error_factors) {
                                    int score = 0;
                                    for (int jj = 0; jj < n; jj++) {
                                        if (M(ei, jj) == 1 && B_best(l, jj) == 1) score++;
                                    }
                                    if (score > best_score) {
                                        best_score = score;
                                        best_factor = l;
                                    }
                                }
                                if (best_factor >= 0) {
                                    A_current(ei, best_factor) = 1;
                                    B_current(best_factor, ej) = 1;
                                }
                            } else {
                                // Erreur: devrait être 0, mais c'est 1
                                // Solution: désactiver TOUS les facteurs pour cette cellule (flip aléatoire)
                                for (int l = 0; l < k; l++) {
                                    if (A_best(ei, l) == 1 && B_best(l, ej) == 1) {
                                        // Flip l'un des deux avec probabilité 50%
                                        if (prob_dist(restart_rng) < 0.5) {
                                            A_current(ei, l) = 0;
                                        } else {
                                            B_current(l, ej) = 0;
                                        }
                                        break;  // Un seul flip par erreur
                                    }
                                }
                            }
                        }
                        
                        if (verbose) {
                            std::cout << "  -> [S8] Correction intelligente de " << errors.size() 
                                      << " erreurs (analyse patterns)" << std::endl;
                        }
                    } else if (strategy == 1) {
                        // Stratégie 1: Perturbation MICRO ciblée (1% sur max 1 facteur - réduit)
                        A_current = A_best;
                        B_current = B_best;
                        
                        int count = 0;
                        for (int l : limited_error_factors) {
                            if (count++ >= 1) break;  // Max 1 facteur (réduit de 2)
                            for (int i = 0; i < m; i++) {
                                if (prob_dist(restart_rng) < 0.01) {  // 1% (réduit de 2%)
                                    A_current(i, l) = 1 - A_current(i, l);
                                }
                            }
                            for (int j = 0; j < n; j++) {
                                if (prob_dist(restart_rng) < 0.01) {  // 1%
                                    B_current(l, j) = 1 - B_current(l, j);
                                }
                            }
                        }
                        
                        if (verbose) {
                            std::cout << "  -> [S1] Perturbation micro 1% sur max 1 facteur" << std::endl;
                        }
                    } else if (strategy == 2) {
                        // Stratégie 2: Perturbation LÉGÈRE ciblée (5% sur max 2 facteurs)
                        A_current = A_best;
                        B_current = B_best;
                        
                        int count = 0;
                        for (int l : limited_error_factors) {
                            if (count++ >= 2) break;  // Max 2 facteurs au lieu de 3
                            for (int i = 0; i < m; i++) {
                                if (prob_dist(restart_rng) < 0.05) {  // 5% au lieu de 10%
                                    A_current(i, l) = 1 - A_current(i, l);
                                }
                            }
                            for (int j = 0; j < n; j++) {
                                if (prob_dist(restart_rng) < 0.05) {  // 5% au lieu de 10%
                                    B_current(l, j) = 1 - B_current(l, j);
                                }
                            }
                        }
                        
                        if (verbose) {
                            std::cout << "  -> [S2] Perturbation légère 5% sur max 2 facteurs" << std::endl;
                        }
                    } else if (strategy == 3) {
                        // Stratégie 3: Flip ciblé UNIQUEMENT sur lignes/colonnes d'erreur
                        // AMÉLIORÉE: 10% au lieu de 15% et flip SEULEMENT 1 cellule par ligne/colonne
                        A_current = A_best;
                        B_current = B_best;
                        
                        // Perturber seulement les cellules directement liées aux erreurs
                        // Flip max 1 cellule par ligne pour réduire l'impact
                        for (int i : error_rows_set) {
                            bool flipped = false;
                            for (int l : limited_error_factors) {
                                if (!flipped && prob_dist(restart_rng) < 0.10) {  // 10% au lieu de 15%
                                    A_current(i, l) = 1 - A_current(i, l);
                                    flipped = true;
                                }
                            }
                        }
                        for (int j : error_cols_set) {
                            bool flipped = false;
                            for (int l : limited_error_factors) {
                                if (!flipped && prob_dist(restart_rng) < 0.10) {  // 10% au lieu de 15%
                                    B_current(l, j) = 1 - B_current(l, j);
                                    flipped = true;
                                }
                            }
                        }
                        
                        if (verbose) {
                            std::cout << "  -> [S3] Flip ciblé 10% (max 1/ligne) sur " << error_rows_set.size() 
                                      << " lignes et " << error_cols_set.size() << " colonnes d'erreur" << std::endl;
                        }
                    } else if (strategy == 4) {
                        // Stratégie 4: Réinitialisation d'UN SEUL facteur le plus problématique
                        A_current = A_best;
                        B_current = B_best;
                        
                        // Prendre le premier facteur problématique seulement
                        if (!limited_error_factors.empty()) {
                            int l = limited_error_factors[0];
                            for (int i = 0; i < m; i++) {
                                A_current(i, l) = dist(restart_rng);
                            }
                            for (int j = 0; j < n; j++) {
                                B_current(l, j) = dist(restart_rng);
                            }
                            if (verbose) {
                                std::cout << "  -> [S4] Réinitialisation du facteur " << l << std::endl;
                            }
                        }
                    } else if (strategy == 5) {
                        // Stratégie 5: Inversion d'UN SEUL facteur problématique
                        A_current = A_best;
                        B_current = B_best;
                        
                        if (!limited_error_factors.empty()) {
                            int l = limited_error_factors[total_full_restarts % limited_error_factors.size()];
                            for (int i = 0; i < m; i++) {
                                A_current(i, l) = 1 - A_current(i, l);
                            }
                            for (int j = 0; j < n; j++) {
                                B_current(l, j) = 1 - B_current(l, j);
                            }
                            if (verbose) {
                                std::cout << "  -> [S5] Inversion du facteur " << l << std::endl;
                            }
                        }
                    } else if (strategy == 6) {
                        // Stratégie 6: Perturbation TRÈS LÉGÈRE des lignes d'erreur (5% - réduit)
                        // AMÉLIORÉE: utilise flip au lieu de random, et seulement 5%
                        A_current = A_best;
                        B_current = B_best;
                        
                        int flips = 0;
                        for (int i : error_rows_set) {
                            for (int l : limited_error_factors) {
                                if (prob_dist(restart_rng) < 0.05) {  // 5% au lieu de 10%
                                    A_current(i, l) = 1 - A_current(i, l);  // flip au lieu de random
                                    flips++;
                                }
                            }
                        }
                        
                        if (verbose) {
                            std::cout << "  -> [S6] Flip 5% sur " << error_rows_set.size() 
                                      << " lignes d'erreur (" << flips << " flips)" << std::endl;
                        }
                    } else if (strategy == 7) {
                        // Stratégie 7: Perturbation TRÈS LÉGÈRE des colonnes d'erreur (10%)
                        A_current = A_best;
                        B_current = B_best;
                        
                        for (int j : error_cols_set) {
                            for (int l : limited_error_factors) {
                                if (prob_dist(restart_rng) < 0.10) {
                                    B_current(l, j) = dist(restart_rng);
                                }
                            }
                        }
                        
                        if (verbose) {
                            std::cout << "  -> [S7] Perturbation 10% des " << error_cols_set.size() 
                                      << " colonnes d'erreur" << std::endl;
                        }
                    } else {
                        // Stratégie 0: Perturbation globale ULTRA-LÉGÈRE (2%)
                        A_current = A_best;
                        B_current = B_best;
                        
                        for (int i = 0; i < m; i++) {
                            for (int l = 0; l < k; l++) {
                                if (prob_dist(restart_rng) < 0.02) {
                                    A_current(i, l) = 1 - A_current(i, l);
                                }
                            }
                        }
                        for (int l = 0; l < k; l++) {
                            for (int j = 0; j < n; j++) {
                                if (prob_dist(restart_rng) < 0.02) {
                                    B_current(l, j) = 1 - B_current(l, j);
                                }
                            }
                        }
                        
                        if (verbose) {
                            std::cout << "  -> [S0] Perturbation ultra-légère globale 2%" << std::endl;
                        }
                    }
                    
                    // === RÉPARATION RAPIDE APRÈS RESTART ===
                    // Compter les erreurs après perturbation
                    Matrix post_computed = A_current.multiply(B_current);
                    int post_restart_errors = 0;
                    for (int ii = 0; ii < m; ii++) {
                        for (int jj = 0; jj < n; jj++) {
                            if (M(ii, jj) >= 0 && post_computed(ii, jj) != M(ii, jj)) {
                                post_restart_errors++;
                            }
                        }
                    }
                    
                    if (verbose) {
                        std::cout << "  [POST-RESTART] " << post_restart_errors << " erreurs" << std::endl;
                    }
                    
                    // Si la perturbation a créé trop d'erreurs, on garde la meilleure solution
                    // Les stratégies sont maintenant assez légères pour que ça n'arrive pas souvent
                    if (post_restart_errors > best_cost * 3) {
                        A_current = A_best;
                        B_current = B_best;
                        if (verbose) {
                            std::cout << "  [ROLLBACK] Trop d'erreurs, retour à best (" << best_cost << ")" << std::endl;
                        }
                    }
                    
                    // Reset de l'escalade - recommencer depuis le début
                    recently_destroyed.clear();
                    iterations_level5 = 0;
                    escalation_level = 0;  // Recommencer depuis le niveau 1!
                    
                    // Activer le mode exploration après restart
                    exploration_mode_iterations = EXPLORATION_MODE_DURATION;
                    adaptive_stagnation_threshold = stagnation_threshold * 2;
                    
                    // Réinitialiser le voisinage de base - MAIS PLUS PETIT pour être rapide
                    rows_to_free = 15;  // Revenir à la taille de base
                    cols_to_free = 15;
                }
                
                stagnation = 0;
            }
        }
        
        result.error_history.push_back(best_cost);
    }
    
    // === 4. FINALISATION ===
    auto end_time = std::chrono::high_resolution_clock::now();
    result.total_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    result.iterations = static_cast<int>(result.iteration_times.size());
    result.final_errors = best_cost;
    result.success = (best_cost == 0);
    result.A_solution = A_best;
    result.B_solution = B_best;
    result.total_stagnations = total_restarts;
    
    if (verbose) {
        std::cout << "\n=== LNS-v3 TERMINÉ ===" << std::endl;
        std::cout << "Erreurs: " << result.initial_errors << " -> " << result.final_errors << std::endl;
        std::cout << "Temps: " << result.total_time << " ms" << std::endl;
        std::cout << "Restarts: " << total_restarts << std::endl;
    }
    
    return result;
}

// === Helpers pour LNS-v3 ===

int BMF::lns_step_fix_A(
    Matrix& A_current, Matrix& B_current,
    const std::vector<int>& cols_to_free
) {
    // A est complètement fixé, on optimise seulement B pour les colonnes sélectionnées
    std::vector<bool> free_cols(n, false);
    for (int c : cols_to_free) if (c >= 0 && c < n) free_cols[c] = true;
    
    // Calculer le coût fixe (colonnes non libérées)
    int fixed_cost = 0;
    Matrix computed = A_current.multiply(B_current);
    
    for (int j = 0; j < n; j++) {
        if (free_cols[j]) continue;  // Colonne libérée, sera optimisée
        for (int i = 0; i < m; i++) {
            if (M(i, j) >= 0 && computed(i, j) != M(i, j)) {
                fixed_cost++;
            }
        }
    }
    
    // Créer solveur pour les colonnes libres
    SATSolver partial_solver;
    VariableManager B_partial(partial_solver);
    
    // Variables pour B[l, j] avec j libre
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            if (free_cols[j]) {
                B_partial.get(l, j);
            }
        }
    }
    
    // Pour chaque cellule (i, j) avec j libre, créer les contraintes
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (!free_cols[j]) continue;
            if (M(i, j) == -1) continue;  // Don't care
            
            // (A×B)[i,j] = OR_l (A[i,l] AND B[l,j])
            // A est fixé, donc on sait quels A[i,l] = 1
            std::vector<int> active_factors;
            for (int l = 0; l < k; l++) {
                if (A_current(i, l) == 1) {
                    active_factors.push_back(l);
                }
            }
            
            if (M(i, j) == 1) {
                // OR des B[l,j] pour l actif doit être 1
                if (active_factors.empty()) {
                    // Impossible de couvrir ce 1 - erreur inévitable
                    fixed_cost++;
                } else {
                    std::vector<int> clause;
                    for (int l : active_factors) {
                        clause.push_back(B_partial.get(l, j));
                    }
                    partial_solver.add_soft_clause(clause, 1);
                }
            } else {
                // M[i,j] = 0: tous les B[l,j] doivent être 0 pour l actif
                for (int l : active_factors) {
                    partial_solver.add_soft_clause({-B_partial.get(l, j)}, 1);
                }
            }
        }
    }
    
    // Résoudre
    if (!partial_solver.solve()) {
        return -1;
    }
    
    int solver_cost = partial_solver.get_cost();
    
    // Extraire solution
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            if (free_cols[j]) {
                auto it = B_partial.vars.find({l, j});
                if (it != B_partial.vars.end()) {
                    B_current(l, j) = partial_solver.getValue(it->second) ? 1 : 0;
                }
            }
        }
    }
    
    return fixed_cost + solver_cost;
}

int BMF::lns_step_fix_B(
    Matrix& A_current, Matrix& B_current,
    const std::vector<int>& rows_to_free
) {
    // B est complètement fixé, on optimise seulement A pour les lignes sélectionnées
    std::vector<bool> free_rows(m, false);
    for (int r : rows_to_free) if (r >= 0 && r < m) free_rows[r] = true;
    
    // Calculer le coût fixe (lignes non libérées)
    int fixed_cost = 0;
    Matrix computed = A_current.multiply(B_current);
    
    for (int i = 0; i < m; i++) {
        if (free_rows[i]) continue;
        for (int j = 0; j < n; j++) {
            if (M(i, j) >= 0 && computed(i, j) != M(i, j)) {
                fixed_cost++;
            }
        }
    }
    
    // Créer solveur pour les lignes libres
    SATSolver partial_solver;
    VariableManager A_partial(partial_solver);
    
    // Variables pour A[i, l] avec i libre
    for (int i = 0; i < m; i++) {
        if (free_rows[i]) {
            for (int l = 0; l < k; l++) {
                A_partial.get(i, l);
            }
        }
    }
    
    // Pour chaque cellule (i, j) avec i libre, créer les contraintes
    for (int i = 0; i < m; i++) {
        if (!free_rows[i]) continue;
        for (int j = 0; j < n; j++) {
            if (M(i, j) == -1) continue;
            
            // (A×B)[i,j] = OR_l (A[i,l] AND B[l,j])
            // B est fixé, donc on sait quels B[l,j] = 1
            std::vector<int> active_factors;
            for (int l = 0; l < k; l++) {
                if (B_current(l, j) == 1) {
                    active_factors.push_back(l);
                }
            }
            
            if (M(i, j) == 1) {
                if (active_factors.empty()) {
                    fixed_cost++;
                } else {
                    std::vector<int> clause;
                    for (int l : active_factors) {
                        clause.push_back(A_partial.get(i, l));
                    }
                    partial_solver.add_soft_clause(clause, 1);
                }
            } else {
                // M[i,j] = 0: tous les A[i,l] doivent être 0 pour l actif
                for (int l : active_factors) {
                    partial_solver.add_soft_clause({-A_partial.get(i, l)}, 1);
                }
            }
        }
    }
    
    // Résoudre
    if (!partial_solver.solve()) {
        return -1;
    }
    
    int solver_cost = partial_solver.get_cost();
    
    // Extraire solution
    for (int i = 0; i < m; i++) {
        if (free_rows[i]) {
            for (int l = 0; l < k; l++) {
                auto it = A_partial.vars.find({i, l});
                if (it != A_partial.vars.end()) {
                    A_current(i, l) = partial_solver.getValue(it->second) ? 1 : 0;
                }
            }
        }
    }
    
    return fixed_cost + solver_cost;
}