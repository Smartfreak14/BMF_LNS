#include "BMFLocalSearch.hpp"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <cmath>
#include <climits>
#include <unordered_set>

BMFLocalSearch::BMFLocalSearch(int k_factors, const Matrix& matrix, unsigned seed)
    : m(matrix.rows), n(matrix.cols), k(k_factors), M(matrix),
      A(m, k, 0), B(k, n, 0),
      // OPTIMISATION: Allocation contigue pour meilleur cache
      count_flat(m * n, 0),
      score_A_flat(m * k, 0),
      score_B_flat(k * n, 0),
      // OPTIMISATION: unordered_set pour O(1) insert/remove
      A_ones_set(k),
      B_ones_set(k),
      rng(seed) {
    // Initialiser les vecteurs legacy pour compatibilité
    count.resize(m, std::vector<int>(n, 0));
    score_A.resize(m, std::vector<int>(k, 0));
    score_B.resize(k, std::vector<int>(n, 0));
    B_ones_cols.resize(k);
    A_ones_rows.resize(k);
}


void BMFLocalSearch::initialize_greedy() {
    // Initialisation gloutonne: pour chaque facteur, couvrir les 1 non couverts
    
    // Matrice des 1 encore non couverts
    std::vector<std::vector<bool>> uncovered(m, std::vector<bool>(n, false));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (M(i, j) == 1) {
                uncovered[i][j] = true;
            }
        }
    }
    
    // Initialiser A et B à 0
    for (int i = 0; i < m; i++)
        for (int l = 0; l < k; l++)
            A(i, l) = 0;
    for (int l = 0; l < k; l++)
        for (int j = 0; j < n; j++)
            B(l, j) = 0;
    
    // Pour chaque facteur
    for (int l = 0; l < k; l++) {
        // Trouver la colonne avec le plus de 1 non couverts
        int best_col = -1;
        int best_count = 0;
        for (int j = 0; j < n; j++) {
            int count = 0;
            for (int i = 0; i < m; i++) {
                if (uncovered[i][j]) count++;
            }
            if (count > best_count) {
                best_count = count;
                best_col = j;
            }
        }
        
        if (best_col == -1 || best_count == 0) break;  // Plus de 1 à couvrir
        
        // B[l, best_col] = 1
        B(l, best_col) = 1;
        
        // Pour chaque ligne i où M[i, best_col] = 1, A[i, l] = 1
        for (int i = 0; i < m; i++) {
            if (M(i, best_col) == 1) {
                A(i, l) = 1;
                // Marquer toutes les cellules couvertes par cette ligne
                for (int j = 0; j < n; j++) {
                    if (B(l, j) == 1 && uncovered[i][j]) {
                        uncovered[i][j] = false;
                    }
                }
            }
        }
        
        // Étendre B[l, :] pour couvrir plus de 1 avec les lignes sélectionnées
        for (int j = 0; j < n; j++) {
            if (B(l, j) == 0) {
                int gain = 0;
                for (int i = 0; i < m; i++) {
                    if (A(i, l) == 1) {
                        if (uncovered[i][j]) gain++;  // 1 non couvert qu'on va couvrir
                        else if (M(i, j) == 0) gain--;  // 0 qu'on va mal couvrir
                    }
                }
                if (gain > 0) {
                    B(l, j) = 1;
                    // Marquer comme couvert
                    for (int i = 0; i < m; i++) {
                        if (A(i, l) == 1 && uncovered[i][j]) {
                            uncovered[i][j] = false;
                        }
                    }
                }
            }
        }
    }
}


void BMFLocalSearch::compute_all_counts() {
    // OPTIMISATION: Utilise la mémoire contigue pour meilleur cache
    // count[i][j] = nombre de l tels que A[i,l] = 1 ET B[l,j] = 1
    std::fill(count_flat.begin(), count_flat.end(), 0);
    
    for (int i = 0; i < m; i++) {
        for (int l = 0; l < k; l++) {
            if (A(i, l) == 1) {
                // Pour chaque j où B[l,j] = 1, incrémenter count[i,j]
                for (int j = 0; j < n; j++) {
                    if (B(l, j) == 1) {
                        count_flat[i * n + j]++;
                    }
                }
            }
        }
    }
    
    // Synchroniser avec la structure legacy
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            count[i][j] = count_flat[i * n + j];
        }
    }
}


void BMFLocalSearch::flip_A(int i, int l) {
    // 1. Mettre à jour count pour toutes les cellules (i, j) où B[l,j] = 1
    int delta = (A(i, l) == 0) ? 1 : -1;
    const int base_idx = i * n; 
    
    for (int j = 0; j < n; j++) {
        if (B(l, j) == 1) {
            count_flat[base_idx + j] += delta;
            count[i][j] += delta;  
        }
    }
    
    // 2. Flipper A[i,l]
    A(i, l) = 1 - A(i, l);
    
    // NOTE: Les étapes 3, 4, 5 de mise à jour des scores non pondérés 
    // sont supprimées car WLS calcule les scores pondérés à la volée.
    // Les scores non pondérés (score_A, score_B) ne sont plus maintenus.
}

void BMFLocalSearch::flip_B(int l, int j) {
    // 1. Mettre à jour count pour toutes les cellules (i, j) où A[i,l] = 1
    int delta = (B(l, j) == 0) ? 1 : -1;
    
    for (int i = 0; i < m; i++) {
        if (A(i, l) == 1) {
            count_flat[i * n + j] += delta;
            count[i][j] += delta;  // Sync legacy
        }
    }
    
    // 2. Flipper B[l,j]
    B(l, j) = 1 - B(l, j);
    
    // NOTE: Les étapes 3, 4, 5 de mise à jour des scores non pondérés 
    // sont supprimées car WLS calcule les scores pondérés à la volée.
    // Les scores non pondérés (score_A, score_B) ne sont plus maintenus.
}

int BMFLocalSearch::count_errors() const {
    int errors = 0;
    for (int i = 0; i < m; i++) {
        const int base_idx = i * n;  // OPTIMISATION: Précalculer l'offset
        for (int j = 0; j < n; j++) {
            if (M(i, j) == -1) continue;
            int computed = (count_flat[base_idx + j] > 0) ? 1 : 0;
            if (computed != M(i, j)) errors++;
        }
    }
    return errors;
}

// ==================== MÉTHODE BASIC====================

 
LocalSearchResult BMFLocalSearch::solve_basic(int /* max_iterations */, bool verbose) {
    LocalSearchResult result;
    result.method_name = "BASIC";
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Initialisation
    compute_all_counts();
    
    int errors = count_errors();
    result.initial_errors = errors;
    result.error_history.push_back(errors);
    
    if (verbose) {
        std::cout << "  [BASIC] Initial errors: " << errors << std::endl;
    }
    
    int best_errors = errors;
    Matrix A_best = A;
    Matrix B_best = B;
    int iter = 0;
    
    // Boucle jusqu'à minimum local (pas de limite d'itérations)
    while (errors > 0) {
        iter++;
        result.iterations = iter;
        
        // Recalculer TOUS les scores à chaque itération (recherche exhaustive)
        // pour garantir de trouver le meilleur flip disponible
        char best_type = ' ';
        int best_idx1 = -1, best_idx2 = -1;
        int best_score = 0;
        
        // Chercher le meilleur flip dans A
        for (int i = 0; i < m; i++) {
            for (int l = 0; l < k; l++) {
                int score = compute_score_A(i, l);
                if (score > best_score) {
                    best_score = score;
                    best_type = 'A';
                    best_idx1 = i;
                    best_idx2 = l;
                }
            }
        }
        
        // Chercher le meilleur flip dans B
        for (int l = 0; l < k; l++) {
            for (int j = 0; j < n; j++) {
                int score = compute_score_B(l, j);
                if (score > best_score) {
                    best_score = score;
                    best_type = 'B';
                    best_idx1 = l;
                    best_idx2 = j;
                }
            }
        }
        
        // Si aucun flip positif, on a atteint un minimum local
        if (best_type == ' ' || best_score <= 0) {
            if (verbose) {
                std::cout << "  [BASIC] Minimum local atteint à iter " << iter 
                          << " avec " << errors << " erreurs" << std::endl;
            }
            result.stop_reason = "local_minimum";
            break;
        }
        
        // Effectuer le flip
        if (best_type == 'A') {
            flip_A(best_idx1, best_idx2);
        } else {
            flip_B(best_idx1, best_idx2);
        }
        
        errors = count_errors();
        result.error_history.push_back(errors);
        
        if (errors < best_errors) {
            best_errors = errors;
            A_best = A;
            B_best = B;
            
            if (verbose) {
                auto now = std::chrono::high_resolution_clock::now();
                double elapsed = std::chrono::duration<double, std::milli>(now - start_time).count();
                std::cout << "  [BASIC] Iter " << iter << ": " << errors 
                          << " erreurs (gain=" << best_score << ", flip=" << best_type 
                          << ", t=" << std::fixed << std::setprecision(0) << elapsed << "ms)" << std::endl;
            }
        }
        
        if (errors == 0) {
            result.success = true;
            break;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.total_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    result.final_errors = best_errors;
    result.A_solution = A_best;
    result.B_solution = B_best;
    result.success = (best_errors == 0);
    
    return result;
}

int BMFLocalSearch::compute_score_A(int i, int l) {
    // Score = gain d'erreurs si on flippe A[i,l]
    // gain > 0 = amélioration (moins d'erreurs après flip)
    int gain = 0;
    const int base_idx = i * n; 
    
    if (A(i, l) == 0) {
        // Flip 0 -> 1
        for (int j = 0; j < n; j++) {
            if (B(l, j) == 1) {
                if (M(i, j) == -1) continue;  
                
                int old_count = count_flat[base_idx + j];
                if (old_count == 0) {
                    if (M(i, j) == 1) gain++;   // Erreur corrigée
                    else gain--;                 // Erreur créée
                }
            }
        }
    } else {
        // Flip 1 -> 0
        for (int j = 0; j < n; j++) {
            if (B(l, j) == 1) {
                if (M(i, j) == -1) continue;
                
                int old_count = count_flat[base_idx + j];
                if (old_count == 1) {
                    if (M(i, j) == 1) gain--;   // Erreur créée
                    else gain++;                 // Erreur corrigée
                }
            }
        }
    }
    
    return gain;
}

int BMFLocalSearch::compute_score_B(int l, int j) {
    int gain = 0;
    
    if (B(l, j) == 0) {
        // Flip 0 -> 1
        for (int i = 0; i < m; i++) {
            if (A(i, l) == 1) {
                if (M(i, j) == -1) continue;
                
                int old_count = count_flat[i * n + j];  // OPTIMISATION: accès contigu
                if (old_count == 0) {
                    if (M(i, j) == 1) gain++;
                    else gain--;
                }
            }
        }
    } else {
        // Flip 1 -> 0
        for (int i = 0; i < m; i++) {
            if (A(i, l) == 1) {
                if (M(i, j) == -1) continue;
                
                int old_count = count_flat[i * n + j];
                if (old_count == 1) {
                    if (M(i, j) == 1) gain--;
                    else gain++;
                }
            }
        }
    }
    
    return gain;
}


std::tuple<char, int, int, int> BMFLocalSearch::find_best_flip() {
    char best_type = ' ';
    int best_idx1 = -1, best_idx2 = -1;
    int best_score = 0;  // On cherche le meilleur score POSITIF
    
    // OPTIMISATION: Utilise la mémoire contiguë pour meilleur cache locality
    // Chercher dans A
    for (int i = 0; i < m; i++) {
        const int base_idx = i * k;
        for (int l = 0; l < k; l++) {
            int score = score_A_flat[base_idx + l];
            if (score > best_score) {
                best_score = score;
                best_type = 'A';
                best_idx1 = i;
                best_idx2 = l;
            }
        }
    }
    
    // Chercher dans B
    for (int l = 0; l < k; l++) {
        const int base_idx = l * n;
        for (int j = 0; j < n; j++) {
            int score = score_B_flat[base_idx + j];
            if (score > best_score) {
                best_score = score;
                best_type = 'B';
                best_idx1 = l;
                best_idx2 = j;
            }
        }
    }
    
    return {best_type, best_idx1, best_idx2, best_score};
}


void BMFLocalSearch::compute_all_scores() {
    for (int i = 0; i < m; i++) {
        for (int l = 0; l < k; l++) {
            int score = compute_score_A(i, l);
            score_A_flat[i * k + l] = score;
            score_A[i][l] = score;  
        }
    }
    
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            int score = compute_score_B(l, j);
            score_B_flat[l * n + j] = score;
            score_B[l][j] = score;  
        }
    }
}
// ==================== INTERFACE LEGACY ====================

LocalSearchResult BMFLocalSearch::solve(int max_iterations, bool verbose) {
    // Par défaut, utiliser BASIC
    return solve_basic(max_iterations, verbose);
}



// ==================== WEIGHTED LOCAL SEARCH ====================

double BMFLocalSearch::compute_weighted_score_A(int i, int l) {
    // Score pondéré = somme des gains * poids des cellules
    double gain = 0.0;
    const int base_idx = i * n;
    
    if (A(i, l) == 0) {
        // Flip 0 -> 1
        for (int j = 0; j < n; j++) {
            if (B(l, j) == 1) {
                if (M(i, j) == -1) continue;
                
                int old_count = count_flat[base_idx + j];
                double w = weights_flat[base_idx + j];
                if (old_count == 0) {
                    if (M(i, j) == 1) gain += w;   // Erreur corrigée (pondérée)
                    else gain -= w;                 // Erreur créée (pondérée)
                }
            }
        }
    } else {
        // Flip 1 -> 0
        for (int j = 0; j < n; j++) {
            if (B(l, j) == 1) {
                if (M(i, j) == -1) continue;
                
                int old_count = count_flat[base_idx + j];
                double w = weights_flat[base_idx + j];
                if (old_count == 1) {
                    if (M(i, j) == 1) gain -= w;   // Erreur créée
                    else gain += w;                 // Erreur corrigée
                }
            }
        }
    }
    
    return gain;
}

double BMFLocalSearch::compute_weighted_score_B(int l, int j) {
    double gain = 0.0;
    
    if (B(l, j) == 0) {
        for (int i = 0; i < m; i++) {
            if (A(i, l) == 1) {
                if (M(i, j) == -1) continue;
                
                int old_count = count_flat[i * n + j];
                double w = weights_flat[i * n + j];
                if (old_count == 0) {
                    if (M(i, j) == 1) gain += w;
                    else gain -= w;
                }
            }
        }
    } else {
        for (int i = 0; i < m; i++) {
            if (A(i, l) == 1) {
                if (M(i, j) == -1) continue;
                
                int old_count = count_flat[i * n + j];
                double w = weights_flat[i * n + j];
                if (old_count == 1) {
                    if (M(i, j) == 1) gain -= w;
                    else gain += w;
                }
            }
        }
    }
    
    return gain;
}

std::tuple<char, int, int, double> BMFLocalSearch::find_best_weighted_flip(int current_errors) {
    char best_type = ' ';
    int best_idx1 = -1, best_idx2 = -1;
    double best_score = 0.0;
    
    // SEUIL ADAPTATIF
    double threshold = std::max(2.0, std::min(5.0, 5000.0 / std::max(1, current_errors)));
    
    // FOCUS ERREURS: 70% guidé par erreurs, 30% random
    int total_samples = 400;  // Légèrement réduit
    int error_guided = total_samples * 7 / 10;
    int random_samples = total_samples - error_guided;
    
    std::uniform_int_distribution<> row_dist(0, m - 1);
    std::uniform_int_distribution<> col_k_dist(0, k - 1);
    std::uniform_int_distribution<> col_n_dist(0, n - 1);
    
    // 1. Échantillons GUIDÉS PAR ERREURS (plus efficace)
    // Chercher des cellules en erreur et tester les flips les affectant
    int error_found = 0;
    for (int attempt = 0; attempt < error_guided * 3 && error_found < error_guided; attempt++) {
        int i = row_dist(rng);
        int j = col_n_dist(rng);
        if (M(i, j) == -1) continue;
        
        int computed = (count_flat[i * n + j] > 0) ? 1 : 0;
        if (computed != M(i, j)) {
            // Cellule en erreur trouvée - tester un flip qui l'affecte
            error_found++;
            int l = col_k_dist(rng);
            
            // Alterner entre A et B
            if (error_found % 2 == 0) {
                double score = compute_weighted_score_A(i, l);
                if (score > best_score) {
                    best_score = score;
                    best_type = 'A';
                    best_idx1 = i;
                    best_idx2 = l;
                    if (score > threshold) return {best_type, best_idx1, best_idx2, best_score};
                }
            } else {
                double score = compute_weighted_score_B(l, j);
                if (score > best_score) {
                    best_score = score;
                    best_type = 'B';
                    best_idx1 = l;
                    best_idx2 = j;
                    if (score > threshold) return {best_type, best_idx1, best_idx2, best_score};
                }
            }
        }
    }
    
    // 2. Échantillons RANDOM (exploration)
    for (int s = 0; s < random_samples; s++) {
        if (s % 2 == 0) {
            int i = row_dist(rng);
            int l = col_k_dist(rng);
            double score = compute_weighted_score_A(i, l);
            if (score > best_score) {
                best_score = score;
                best_type = 'A';
                best_idx1 = i;
                best_idx2 = l;
                if (score > threshold) return {best_type, best_idx1, best_idx2, best_score};
            }
        } else {
            int l = col_k_dist(rng);
            int j = col_n_dist(rng);
            double score = compute_weighted_score_B(l, j);
            if (score > best_score) {
                best_score = score;
                best_type = 'B';
                best_idx1 = l;
                best_idx2 = j;
                if (score > threshold) return {best_type, best_idx1, best_idx2, best_score};
            }
        }
    }
    
    return {best_type, best_idx1, best_idx2, best_score};
}

LocalSearchResult BMFLocalSearch::solve_weighted(int max_iterations, 
                                                  double penalty_increment,
                                                  int max_stagnation,
                                                  bool verbose,
                                                  std::atomic<bool>* stop_flag) {
    LocalSearchResult result;
    result.method_name = "WEIGHTED";
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Lambda pour vérifier l'interruption
    auto is_interrupted = [&stop_flag]() -> bool {
        return stop_flag && stop_flag->load();
    };
    
    // Initialiser les poids à 1.0
    weights_flat.assign(m * n, 1.0);
    
    // Initialisation des structures de comptage
    compute_all_counts();
    
    int errors = count_errors();
    result.initial_errors = errors;
    result.error_history.push_back(errors);
    
    if (verbose) {
        std::cout << "  [WEIGHTED] Initial errors: " << errors << std::endl;
    }
    
    int best_errors = errors;
    Matrix A_best = A;
    Matrix B_best = B;
    
    int stagnation = 0;
    int total_penalty_rounds = 0;
    int rounds_without_improvement = 0;
    int iter = 0;
    int restarts = 0;
    const int MAX_RESTARTS = 2;  // Nombre max de redémarrages
    const int RESTART_THRESHOLD = 100;  // Rounds sans amélioration avant restart
    
    while (iter < max_iterations && errors > 0 && !is_interrupted()) {
        // Phase 1: Recherche locale classique jusqu'à stagnation
        auto [flip_type, idx1, idx2, score] = find_best_weighted_flip(errors);
        
        if (flip_type == ' ' || score <= 1e-9) {
            // Minimum local atteint - appliquer les pénalités
            stagnation++;
            
            if (stagnation >= max_stagnation) {
                // Augmenter les poids des cellules en erreur
                total_penalty_rounds++;
                int penalized = 0;
                
                // Augmenter plus agressivement après plusieurs rounds
                double actual_penalty = penalty_increment * (3.0 + 0.1 * (total_penalty_rounds / 20));
                
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < n; j++) {
                        if (M(i, j) == -1) continue;
                        
                        int computed = (count_flat[i * n + j] > 0) ? 1 : 0;
                        if (computed != M(i, j)) {
                            // Cette cellule est en erreur - augmenter son poids
                            weights_flat[i * n + j] += actual_penalty;
                            penalized++;
                        }
                    }
                }
                
                if (verbose && total_penalty_rounds % 10 == 1) {
                    std::cout << "  [WEIGHTED] Penalty round " << total_penalty_rounds 
                              << ": " << penalized << " cells penalized (penalty=" 
                              << std::fixed << std::setprecision(2) << actual_penalty 
                              << "), best=" << best_errors << std::endl;
                }
                
                stagnation = 0;
                rounds_without_improvement++;
                
                // RESTART: Si stagnation trop longue, réinitialiser les poids et perturber
                if (rounds_without_improvement >= RESTART_THRESHOLD && restarts < MAX_RESTARTS) {
                    restarts++;
                    if (verbose) {
                        std::cout << "  [WEIGHTED] RESTART " << restarts << "/" << MAX_RESTARTS 
                                  << ": perturbation guidée par les erreurs" << std::endl;
                    }
                    
                    // Réinitialiser les poids à 1.0
                    std::fill(weights_flat.begin(), weights_flat.end(), 1.0);
                    
                    // === PERTURBATION GUIDÉE PAR LES ERREURS ===
                    // Au lieu de flipper au hasard, on cible les cellules en erreur
                    
                    // 1. Collecter les cellules en erreur
                    std::vector<std::pair<int, int>> error_cells;
                    for (int i = 0; i < m; i++) {
                        for (int j = 0; j < n; j++) {
                            if (M(i, j) == -1) continue;
                            int computed = (count_flat[i * n + j] > 0) ? 1 : 0;
                            if (computed != M(i, j)) {
                                error_cells.push_back({i, j});
                            }
                        }
                    }
                    
                    if (!error_cells.empty()) {
                        std::uniform_int_distribution<> dist_error(0, error_cells.size() - 1);
                        std::uniform_int_distribution<> dist_k(0, k - 1);
                        std::uniform_int_distribution<> dist_choice(0, 1);
                        
                        // Nombre de perturbations adapté au nombre d'erreurs
                        // Plus d'erreurs = plus de perturbations
                        int perturb_count = std::max(1, std::min((int)error_cells.size() / 5, 20));
                        
                        for (int p = 0; p < perturb_count; p++) {
                            // Choisir une cellule EN ERREUR aléatoirement
                            auto [ei, ej] = error_cells[dist_error(rng)];
                            int l = dist_k(rng);
                            
                            // Flipper dans A OU B (pas les deux!) pour éviter le chaos
                            if (dist_choice(rng) == 0) {
                                flip_A(ei, l);
                            } else {
                                flip_B(l, ej);
                            }
                        }
                    }
                    
                    errors = count_errors();
                    rounds_without_improvement = 0;
                    total_penalty_rounds = 0;
                }
                
                // Arrêter si tous les restarts sont épuisés et toujours pas d'amélioration
                if (rounds_without_improvement >= RESTART_THRESHOLD && restarts >= MAX_RESTARTS) {
                    if (verbose) {
                        std::cout << "  [WEIGHTED] Stopping: exhausted all " << MAX_RESTARTS << " restarts" << std::endl;
                    }
                    break;
                }
                
                // Vérifier l'interruption
                if (is_interrupted()) {
                    if (verbose) {
                        std::cout << "  [WEIGHTED] Stopping: interrupted (Ctrl+C)" << std::endl;
                    }
                    break;
                }
            }
        } else {
            // Effectuer le flip
            if (flip_type == 'A') {
                flip_A(idx1, idx2);
            } else {
                flip_B(idx1, idx2);
            }
            
            errors = count_errors();
            result.error_history.push_back(errors);
            stagnation = 0;
            
            if (errors < best_errors) {
                best_errors = errors;
                A_best = A;
                B_best = B;
                rounds_without_improvement = 0;
                
                if (verbose && (iter % 1000 == 0 || errors == 0)) {
                    std::cout << "  [WEIGHTED] Iter " << iter << ": " << errors 
                              << " errors (penalty rounds: " << total_penalty_rounds << ")" << std::endl;
                }
            }
        }
        
        iter++;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.total_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    result.final_errors = best_errors;
    result.iterations = iter;
    result.A_solution = A_best;
    result.B_solution = B_best;
    result.success = (best_errors == 0);
    result.stop_reason = (best_errors == 0) ? "zero_errors" :
                         is_interrupted() ? "interrupted" :
                         (restarts >= MAX_RESTARTS) ? "max_restarts" : "max_iterations";
    
    if (verbose) {
        std::string stop_emoji = (result.stop_reason == "zero_errors") ? "OK" :
                                 (result.stop_reason == "interrupted") ? "INT" : "REST";
        std::cout << "  [WEIGHTED] Final: " << best_errors << " errors in " 
                  << iter << " iterations, " << total_penalty_rounds << " penalty rounds, "
                  << std::fixed << std::setprecision(1) << result.total_time << " ms"
                  << " (" << stop_emoji << " " << result.stop_reason << ")" << std::endl;
    }
    
    return result;
}

