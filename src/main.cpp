#include "BMF.hpp"
#include "BMFLocalSearch.hpp"
#include "CSVMatrixLoader.hpp"
#include <iostream>
#include <iomanip>
#include <random>
#include <map>
#include <set>
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <array>
#include <csignal>
#include <atomic>
#include <sstream>


const std::string CSV_BASE_DIR = "../data/Real/SimplUni/";
const std::map<std::string, int> CSV_FILES_WITH_K = {
    {"zoo.csv", 25},
    {"wine.csv", 178},
    {"votes.csv", 17},
    {"tumor.csv", 44},
    {"tictactoe.csv", 28},
    {"ThoraricSurgery.csv", 304},
    {"soybean", 90},
    {"studentPerformance.csv", 176},
    {"phishing.csv", 26},
    {"paleo.csv", 139},
    {"nursery.csv", 30},
    {"mushroom.csv", 51},
    {"lymphography.csv", 54},
    {"hepatits.csv", 154},
    {"heart.csv", 270},
    {"iris.csv", 121},
    {"flare.csv", 42},
    {"firewall.csv", 65},
    {"dna.csv", 367},
    {"dermatology.csv", 194},
    {"customer.csv", 277},
    {"cmc.csv", 71},
    {"chess_krvskp.csv", 38},
    {"car.csv", 25},
    {"breast.csv", 90},
    {"balance.csv", 23},
    {"autism.csv", 154},
    {"audio.csv", 200},
    {"apj.csv", 454},
    {"americas_small.csv", 187},
    {"advertisement.csv", 749}
};


const std::string RESULTS_OUTPUT_FILE = "../data/resultats_batch.txt";


struct CSVTestResult {
    std::string filename;
    int k;
    int init_errors;
    int after_ls_errors;
    int final_errors;
    int iterations;
    double ls_time;
    double lns_time;
    double total_time;
    bool success;  
    
    std::string to_string() const {
        std::ostringstream oss;
        oss << filename << " k=" << k << " :\n";
        oss << "=== RÉSULTATS LS + WLS + LNS-v3 ===\n";
        oss << "Erreurs: " << init_errors << " -> " << after_ls_errors 
            << " (LS+WLS) -> " << final_errors << " (LNS-v3)\n";
        oss << "Itérations LNS-v3: " << iterations << "\n";
        oss << std::fixed << std::setprecision(1);
        oss << "Temps LS+WLS: " << ls_time << " ms\n";
        oss << "Temps LNS-v3: " << lns_time << " ms\n";
        oss << "Temps total: " << total_time << " ms\n";
        return oss.str();
    }
};

// ==================== GESTION DES SIGNAUX (Ctrl+C) ====================


std::atomic<bool> g_interrupted(false);


void signal_handler(int signum) {
    if (signum == SIGINT) {
        std::cout << "\n\nINTERRUPTION DÉTECTÉE (Ctrl+C) - Arrêt en cours..." << std::endl;
        g_interrupted.store(true);
    }
}


bool is_interrupted() {
    return g_interrupted.load();
}


void reset_interrupt() {
    g_interrupted.store(false);
}


void save_matrix_csv(const Matrix& M, const std::string& filepath) {
    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir " << filepath << std::endl;
        return;
    }
    
    for (int i = 0; i < M.rows; i++) {
        for (int j = 0; j < M.cols; j++) {
            file << M(i, j);
            if (j < M.cols - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
}


int verify_with_python(const Matrix& M, const Matrix& A, const Matrix& B, 
                       const std::string& method_name, const std::string& csv_filename = "") {
    // Créer un dossier temporaire pour les fichiers
    std::string temp_dir = "/tmp/bmf_verify";
    std::filesystem::create_directories(temp_dir);
    
    // Générer des noms de fichiers uniques
    std::string base_name = csv_filename.empty() ? "matrix" : 
                            std::filesystem::path(csv_filename).stem().string();
    std::string suffix = "_" + method_name;
    // Remplacer les espaces par des underscores
    std::replace(suffix.begin(), suffix.end(), ' ', '_');
    std::replace(suffix.begin(), suffix.end(), '+', '_');
    
    std::string A_path = temp_dir + "/" + base_name + suffix + "_A.csv";
    std::string B_path = temp_dir + "/" + base_name + suffix + "_B.csv";
    std::string M_path = temp_dir + "/" + base_name + "_M.csv";
    
    // Sauvegarder les matrices
    save_matrix_csv(A, A_path);
    save_matrix_csv(B, B_path);
    save_matrix_csv(M, M_path);
    
    // Construire la commande Python
    // Le script verif.py est dans ../src/verif.py par rapport au build
    std::string python_cmd = "python3 ../src/verif.py \"" + A_path + "\" \"" + B_path + "\" \"" + M_path + "\" 2>&1";
    
    // Exécuter et capturer la sortie
    std::array<char, 256> buffer;
    std::string output;
    FILE* pipe = popen(python_cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "  [VERIF] Erreur: Impossible d'exécuter verif.py" << std::endl;
        return -1;
    }
    
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        output += buffer.data();
    }
    
    pclose(pipe);
    
    // Parser la sortie pour extraire l'erreur
    int python_errors = -1;
    std::string error_prefix = "Reconstruction error =";
    size_t pos = output.find(error_prefix);
    if (pos != std::string::npos) {
        try {
            python_errors = std::stoi(output.substr(pos + error_prefix.length()));
        } catch (...) {
            std::cerr << "  [VERIF] Erreur: Impossible de parser la sortie" << std::endl;
        }
    }
    
    // Afficher le résultat de vérification
    if (python_errors >= 0) {
        std::cout << "  [VERIF Python] " << method_name << ": " << python_errors << " erreurs";
        if (python_errors == 0) {
            std::cout << " ";
        }
        std::cout << std::endl;
    } else {
        std::cerr << "  [VERIF] Erreur avec verif.py: " << output << std::endl;
    }
    
    return python_errors;
}

// ==================== GÉNÉRATEUR DE DONNÉES ====================

class DataGenerator {
public:
    std::mt19937 rng;
    
    DataGenerator(unsigned seed = 42) : rng(seed) {}
    
    Matrix random_matrix(int rows, int cols, double density = 0.1) {
        Matrix M(rows, cols, 0);
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (dis(rng) < density) {
                    M(i, j) = 1;
                }
            }
        }
        return M;
    }
    
    std::pair<Matrix, Matrix> generate_rbac_config(int m, int n_objects, int k, double density = 0.3) {
        Matrix A = random_matrix(m, k, density);
        Matrix B = random_matrix(k, n_objects * 2, density);
        return {A, B};
    }
    
    
    std::tuple<Matrix, Matrix, Matrix> generate_bmf_instance(
        int m, int n, int k,
        double target_density = 0.1,
        double missing_rate = 0.0
    ) {
        
        double p = std::sqrt(1.0 - std::pow(1.0 - target_density, 1.0 / k));
        
        std::cout << "  [BMF Instance] m=" << m << ", n=" << n << ", k=" << k << std::endl;
        std::cout << "  Target density: " << (target_density * 100) << "%" << std::endl;
        std::cout << "  Entry probability p: " << p << std::endl;
        if (missing_rate > 0) {
            std::cout << "  Missing rate: " << (missing_rate * 100) << "%" << std::endl;
        }
        
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        // 1. Générer A (m×k) avec probabilité p
        Matrix A_true(m, k, 0);
        for (int i = 0; i < m; i++) {
            for (int l = 0; l < k; l++) {
                if (dis(rng) < p) {
                    A_true(i, l) = 1;
                }
            }
        }
        
        // 2. Générer B (k×n) avec probabilité p
        Matrix B_true(k, n, 0);
        for (int l = 0; l < k; l++) {
            for (int j = 0; j < n; j++) {
                if (dis(rng) < p) {
                    B_true(l, j) = 1;
                }
            }
        }
        
        // 3. Calculer M = A x B 
        Matrix M = A_true.multiply(B_true);
        
        int ones_count = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (M(i, j) == 1) ones_count++;
            }
        }
        double actual_density = static_cast<double>(ones_count) / (m * n);
        std::cout << "  Actual density: " << (actual_density * 100) << "% (" << ones_count << "/" << (m*n) << " ones)" << std::endl;
        
        // 4. retirer des entrées et les mettre à -1
        if (missing_rate > 0) {
            int missing_count = 0;
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    if (dis(rng) < missing_rate) {
                        M(i, j) = -1;
                        missing_count++;
                    }
                }
            }
            std::cout << "  Removed entries: " << missing_count << " (" 
                      << (100.0 * missing_count / (m * n)) << "%)" << std::endl;
        }
        
        return {M, A_true, B_true};
    }
    
    // generer seulement M sans retourner A et B 
    Matrix generate_bmf_matrix(
        int m, int n, int k,
        double target_density = 0.1,
        double missing_rate = 0.0
    ) {
        auto [M, A_true, B_true] = generate_bmf_instance(m, n, k, target_density, missing_rate);
        return M;
    }
    
    // statistique d'une matice M
    static void print_matrix_stats(const Matrix& M, const std::string& name = "M") {
        int total = M.rows * M.cols;
        int count_zeros = 0;
        int count_ones = 0;
        int count_missing = 0;
        
        for (int i = 0; i < M.rows; i++) {
            for (int j = 0; j < M.cols; j++) {
                if (M(i, j) == 0) count_zeros++;
                else if (M(i, j) == 1) count_ones++;
                else if (M(i, j) == -1) count_missing++;
            }
        }
        
        std::cout << "=== Statistiques de " << name << " (" << M.rows << "x" << M.cols << ") ===" << std::endl;
        std::cout << "  Total: " << total << " cellules" << std::endl;
        std::cout << "  Ones (1):    " << std::setw(6) << count_ones << " (" 
                  << std::fixed << std::setprecision(2) << (100.0 * count_ones / total) << "%)" << std::endl;
        std::cout << "  Zeros (0):   " << std::setw(6) << count_zeros << " (" 
                  << std::fixed << std::setprecision(2) << (100.0 * count_zeros / total) << "%)" << std::endl;
        std::cout << "  Missing (-1):" << std::setw(6) << count_missing << " (" 
                  << std::fixed << std::setprecision(2) << (100.0 * count_missing / total) << "%)" << std::endl;
        
        // Densité effective (excluant les -1)
        if (count_missing > 0) {
            int known = count_ones + count_zeros;
            std::cout << "  Densité effective (1 parmi connus): " 
                      << std::fixed << std::setprecision(2) << (100.0 * count_ones / known) << "%" << std::endl;
        }
    }
};



//Test la factorisation BMF sur un fichier CSV avec GREEDY + WLS + LNS-v3
CSVTestResult test_csv_partial(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    // Installer le gestionnaire de signal pour Ctrl+C
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== COMPARAISON RANDOM vs LS + PARTIAL + HYBRID (CSV) ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::string filename = csv_file;
    
    // Si c'est juste un nom de fichier, ajouter le chemin data/
    if (filename.find('/') == std::string::npos) {
        filename = "../data/" + filename;
        if (filename.find(".csv") == std::string::npos) {
            filename += ".csv";
        }
    }
    
    std::cout << "Chargement de: " << filename << std::endl;
    
    Matrix M = CSVMatrixLoader::loadFromCSV(filename);
    DataGenerator::print_matrix_stats(M, "Matrice CSV");
    if (M.rows == 0 || M.cols == 0) {
        std::cerr << "Erreur: impossible de charger la matrice" << std::endl;
        return result;  // Retourner résultat vide avec success=false
    }
    
    int m = M.rows;
    int n = M.cols;
    
    // Déterminer k automatiquement si non spécifié
    int k = k_value;
    if (k <= 0) {
        k = std::max(5, std::min(50, static_cast<int>(std::log2(std::min(m, n)) * 2)));
        std::cout << "  k automatique: " << k << " (basé sur dimensions " << m << "x" << n << ")" << std::endl;
    } else {
        std::cout << "  k spécifié: " << k << std::endl;
    }
    
    // Afficher le timeout si configuré
    if (timeout_seconds > 0) {
        std::cout << "  Timeout LNS-v3: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "   Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    
    // ==================== TEST : LS + WLS + LNS-v3 (Alternance A/B + Anti-stagnation) ====================
    std::cout << "\n" << std::string(60, '-') << std::endl;
    std::cout << "--- GREEDY + WLS + LNS-v3 ---" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    int v3_init_errors = 0;
    int v3_after_ls_errors = 0;
    int v3_final_errors = 0;
    double v3_ls_time = 0;
    double v3_lns_time = 0;
    double v3_total_time = 0;
    int v3_iterations = 0;
    
    // Stocker les solutions finales pour vérification
    Matrix v3_A_final(0, 0), v3_B_final(0, 0);
    
    {
        auto total_start = std::chrono::high_resolution_clock::now();
        
        // Phase 1: Initialisation GREEDY seulement (pas de BASIC)
        std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
        auto ls_start = std::chrono::high_resolution_clock::now();
        
        BMFLocalSearch ls(k, M, 42);  // seed harmonisé avec --alt
        ls.initialize_greedy();
        
        // Calculer les erreurs après greedy
        ls.compute_all_counts();
        v3_init_errors = ls.count_errors();
        
        std::cout << "  [GREEDY] Erreurs initiales: " << v3_init_errors << std::endl;
        
        // Phase 1b: Passer DIRECTEMENT à WEIGHTED (sans BASIC)
        LocalSearchResult ls_result;
        ls_result.initial_errors = v3_init_errors;
        ls_result.final_errors = v3_init_errors;
        ls_result.A_solution = ls.A;
        ls_result.B_solution = ls.B;
        
        if (v3_init_errors > 0) {
            std::cout << "\n[Phase 1a] Local Search (directement après GREEDY)..." << std::endl;
            
            // Appliquer WEIGHTED directement sur la solution greedy
            LocalSearchResult weighted_result = ls.solve(
                5000,   // max_iterations
                true       // verbose
            );
            
            ls_result = weighted_result;
        }
        
        auto ls_end = std::chrono::high_resolution_clock::now();
        v3_ls_time = std::chrono::duration<double, std::milli>(ls_end - ls_start).count();
        v3_after_ls_errors = ls_result.final_errors;
        
        std::cout << "  GREEDY+WLS Final: " << v3_init_errors << " -> " << v3_after_ls_errors 
                  << " erreurs en " << std::fixed << std::setprecision(1) << v3_ls_time << " ms" 
                  << " (iterations: " << ls_result.iterations << ")" << std::endl;
        
        // Phase 2: LNS-v3 avec Alternance A/B + Anti-stagnation
        if (v3_after_ls_errors > 0) {
            std::cout << "\n[Phase 2] LNS-v3 (Alternance A/B + Anti-stagnation)..." << std::endl;
            
            auto lns_start = std::chrono::high_resolution_clock::now();
            
            BMF bmf(k, M);
            
            // Récupérer la solution du LS
            bmf.set_initial_solution(ls_result.A_solution, ls_result.B_solution);
            
            // Paramètres LNS-v3
            double neighborhood_ratio = 0.05;  // 5% des dimensions
            int max_iter_v3 = 10000;          // Beaucoup d'itérations (timeout contrôle l'arrêt)
            int stagnation_threshold = 15;
            
            // Convertir le timeout en millisecondes
            double timeout_ms = timeout_seconds * 1000.0;
            
            auto lns_result = bmf.solve_lns_v3(
                max_iter_v3,          // max_iterations
                neighborhood_ratio,   // neighborhood_size
                42,                   // seed
                stagnation_threshold, // stagnation_threshold
                true,                 // verbose
                timeout_ms,           // timeout_ms
                &g_interrupted        // stop_flag (pour Ctrl+C)
            );
            
            auto lns_end = std::chrono::high_resolution_clock::now();
            v3_lns_time = std::chrono::duration<double, std::milli>(lns_end - lns_start).count();
            v3_final_errors = lns_result.final_errors;
            v3_iterations = lns_result.iterations;
            
            // Stocker les solutions finales de LNS-v3
            v3_A_final = lns_result.A_solution;
            v3_B_final = lns_result.B_solution;
            
            // Afficher la raison de l'arrêt
            std::string stop_emoji = (lns_result.stop_reason == "zero_errors") ? "OK" :
                                     (lns_result.stop_reason == "timeout") ? "TIMEOUT" :
                                     (lns_result.stop_reason == "interrupted") ? "INT" : "REST";
            
            std::cout << "  LNS-v3: " << v3_after_ls_errors << " -> " << v3_final_errors 
                      << " erreurs en " << v3_iterations << " iters (" 
                      << std::fixed << std::setprecision(1) << v3_lns_time << " ms)" << std::endl;
            std::cout << "  Arrêt: " << stop_emoji << " " << lns_result.stop_reason << std::endl;
        } else {
            v3_final_errors = 0;
            v3_iterations = 0;
            
            // Stocker les solutions finales de WLS (pas de LNS)
            v3_A_final = ls_result.A_solution;
            v3_B_final = ls_result.B_solution;
            
            std::cout << "\n[Phase 2] Pas besoin de LNS-v3, déjà 0 erreurs!" << std::endl;
        }
        
        auto total_end = std::chrono::high_resolution_clock::now();
        v3_total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
        
        std::cout << "\n=== RÉSULTATS LS + WLS + LNS-v3 ===" << std::endl;
        std::cout << "Erreurs: " << v3_init_errors << " -> " << v3_after_ls_errors 
                  << " (LS+WLS) -> " << v3_final_errors << " (LNS-v3)" << std::endl;
        std::cout << "Itérations LNS-v3: " << v3_iterations << std::endl;
        std::cout << "Temps LS+WLS: " << std::fixed << std::setprecision(1) << v3_ls_time << " ms" << std::endl;
        std::cout << "Temps LNS-v3: " << std::fixed << std::setprecision(1) << v3_lns_time << " ms" << std::endl;
        std::cout << "Temps total: " << std::fixed << std::setprecision(1) << v3_total_time << " ms" << std::endl;
        
        // Remplir le résultat
        result.k = k;  // Mettre à jour avec la valeur réelle utilisée
        result.init_errors = v3_init_errors;
        result.after_ls_errors = v3_after_ls_errors;
        result.final_errors = v3_final_errors;
        result.iterations = v3_iterations;
        result.ls_time = v3_ls_time;
        result.lns_time = v3_lns_time;
        result.total_time = v3_total_time;
        result.success = true;
    }
    
    // ==================== COMPARAISON FINALE ====================
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== COMPARAISON FINALE ===" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "Fichier: " << CSVMatrixLoader::getBaseName(filename) << std::endl;
    std::cout << "Dimensions: " << m << "x" << n << ", k=" << k << std::endl;
    std::cout << std::endl;
    
    std::cout << std::setw(30) << "Méthode" 
              << std::setw(15) << "Init" 
              << std::setw(15) << "Après WLS"
              << std::setw(15) << "Final" 
              << std::setw(15) << "Temps (ms)" << std::endl;
    std::cout << std::string(90, '-') << std::endl;
    
    std::cout << std::setw(30) << "GREEDY + WLS + LNS-v3" 
              << std::setw(15) << v3_init_errors 
              << std::setw(15) << v3_after_ls_errors
              << std::setw(15) << v3_final_errors
              << std::setw(15) << std::fixed << std::setprecision(1) << v3_total_time << std::endl;
    
    std::cout << std::string(90, '-') << std::endl;
    
    std::cout << "\n RÉSULTAT: " << v3_final_errors << " erreurs finales" << std::endl;
    std::cout << "   Temps total: " << std::fixed << std::setprecision(1) << v3_total_time << " ms" << std::endl;
    
    // ==================== VÉRIFICATION PYTHON ====================
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== VÉRIFICATION PYTHON (verif.py) ===" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    int python_errors = verify_with_python(M, v3_A_final, v3_B_final, "GREEDY_WLS_LNS-v3", filename);
    
    if (python_errors >= 0) {
        if (python_errors == v3_final_errors) {
            std::cout << " Vérification OK: C++ (" << v3_final_errors << ") == Python (" << python_errors << ")" << std::endl;
        } else {
            std::cout << " Différence: C++ (" << v3_final_errors << ") != Python (" << python_errors << ")" << std::endl;
        }
    }
    
    return result;
}


// ==================== TEST ALTERNANCE WLS <-> LNS-v3 ====================

/**
 * @brief Test la factorisation BMF avec alternance multiple WLS <-> LNS-v3
 * 
 * @param csv_file Chemin vers le fichier CSV
 * @param k_value Nombre de facteurs (0 = automatique)
 * @param timeout_seconds Timeout en secondes (0 = pas de timeout)
 * @return CSVTestResult contenant les résultats du test
 */
CSVTestResult test_csv_alternating(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    // Installer le gestionnaire de signal pour Ctrl+C
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== GREEDY + (LS <-> LNS-v3)* ALTERNANCE ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::string filename = csv_file;
    
    // Si c'est juste un nom de fichier, ajouter le chemin data/
    if (filename.find('/') == std::string::npos) {
        filename = "../data/" + filename;
        if (filename.find(".csv") == std::string::npos) {
            filename += ".csv";
        }
    }
    
    std::cout << "Chargement de: " << filename << std::endl;
    
    Matrix M = CSVMatrixLoader::loadFromCSV(filename);
    DataGenerator::print_matrix_stats(M, "Matrice CSV");
    if (M.rows == 0 || M.cols == 0) {
        std::cerr << "Erreur: impossible de charger la matrice" << std::endl;
        return result;
    }
    
    int m = M.rows;
    int n = M.cols;
    
    // Déterminer k automatiquement si non spécifié
    int k = k_value;
    if (k <= 0) {
        k = std::max(5, std::min(50, static_cast<int>(std::log2(std::min(m, n)) * 2)));
        std::cout << "  k automatique: " << k << " (basé sur dimensions " << m << "x" << n << ")" << std::endl;
    } else {
        std::cout << "  k spécifié: " << k << std::endl;
    }
    
    // Afficher le timeout si configuré
    if (timeout_seconds > 0) {
        std::cout << "  Timeout: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "  Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    // ==================== ALTERNANCE WLS <-> LNS-v3 ====================
    std::cout << "\n" << std::string(60, '-') << std::endl;
    std::cout << "--- GREEDY + (LS <-> LNS-v3)* Alternance ---" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    int v3_init_errors = 0;
    int v3_final_errors = 0;
    double v3_ls_time = 0;
    double v3_lns_time = 0;
    double v3_total_time = 0;
    int v3_iterations = 0;
    int total_alternations = 0;
    
    // Stocker les solutions finales pour vérification
    Matrix v3_A_final(0, 0), v3_B_final(0, 0);
    
    {
        auto total_start = std::chrono::high_resolution_clock::now();
        
        // Convertir le timeout en millisecondes
        double timeout_ms = timeout_seconds * 1000.0;
        auto deadline = (timeout_ms > 0) ? 
            std::chrono::high_resolution_clock::now() + std::chrono::milliseconds(static_cast<long long>(timeout_ms)) :
            std::chrono::high_resolution_clock::time_point::max();
        
        // Lambda pour vérifier le temps restant
        auto time_remaining_ms = [&]() -> double {
            if (timeout_ms <= 0) return std::numeric_limits<double>::max();
            auto now = std::chrono::high_resolution_clock::now();
            return std::chrono::duration<double, std::milli>(deadline - now).count();
        };
        
        // Phase 1: Initialisation GREEDY
        std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
        
        BMFLocalSearch ls(k, M, 42);
        ls.initialize_greedy();
        
        // Calculer les erreurs après greedy
        ls.compute_all_counts();
        v3_init_errors = ls.count_errors();
        
        std::cout << "  [GREEDY] Erreurs initiales: " << v3_init_errors << std::endl;
        
        // Initialiser les solutions courantes
        Matrix A_current = ls.A;
        Matrix B_current = ls.B;
        int current_errors = v3_init_errors;
        int best_errors = current_errors;
        Matrix A_best = A_current;
        Matrix B_best = B_current;
        
        // === BOUCLE D'ALTERNANCE WLS <-> LNS-v3 ===
        const int MAX_ALTERNATIONS = 10;  // Maximum d'alternances
        const int MAX_STAGNATION_ALTERNATIONS = 5;  // Arrêter si pas d'amélioration après N alternances
        int stagnation_alternations = 0;
        int errors_before_alternation = current_errors;
        
        while (current_errors > 0 && 
               total_alternations < MAX_ALTERNATIONS && 
               !is_interrupted() &&
               time_remaining_ms() > 1000) {  // Au moins 1 seconde restante
            
            total_alternations++;
            errors_before_alternation = current_errors;
            
            std::cout << "\n" << std::string(50, '=') << std::endl;
            std::cout << "=== ALTERNANCE " << total_alternations << "/" << MAX_ALTERNATIONS 
                      << " (erreurs: " << current_errors << ") ===" << std::endl;
            std::cout << std::string(50, '=') << std::endl;
            
            // === PHASE WLS ===
            std::cout << "\n[LS] Local Search..." << std::endl;
            auto wls_start = std::chrono::high_resolution_clock::now();
            
            // Créer un nouveau LocalSearch avec la solution courante
            BMFLocalSearch ls_wls(k, M, 42 + total_alternations);
            ls_wls.A = A_current;
            ls_wls.B = B_current;
            ls_wls.compute_all_counts();
            // NOTE: compute_all_scores() supprimé - WLS calcule les scores pondérés à la volée
            
            LocalSearchResult wls_result = ls_wls.solve(
                5000,    // max_iterations (beaucoup)
                true       // verbose (moins verbeux dans la boucle)
            );
            
            auto wls_end = std::chrono::high_resolution_clock::now();
            double wls_time = std::chrono::duration<double, std::milli>(wls_end - wls_start).count();
            v3_ls_time += wls_time;
            
            A_current = wls_result.A_solution;
            B_current = wls_result.B_solution;
            current_errors = wls_result.final_errors;
            
            std::cout << "  [LS] " << errors_before_alternation << " -> " << current_errors 
                      << " erreurs (" << std::fixed << std::setprecision(1) << wls_time << " ms)" << std::endl;
            
            if (current_errors < best_errors) {
                best_errors = current_errors;
                A_best = A_current;
                B_best = B_current;
                stagnation_alternations = 0;
            }
            
            if (current_errors == 0) {
                std::cout << "  Factorisation exacte trouvée par LS!" << std::endl;
                break;
            }
            
            // Vérifier le temps restant
            if (time_remaining_ms() < 1000) {
                std::cout << "  Timeout proche, arrêt." << std::endl;
                break;
            }
            
            // === PHASE LNS-v3 (avec alternance Fix_A, Fix_B, Joint) ===
            std::cout << "\n[LNS-v3] Large Neighborhood Search (Fix_A/Fix_B/Joint)..." << std::endl;
            auto lns_start = std::chrono::high_resolution_clock::now();
            
            BMF bmf(k, M);
            bmf.set_initial_solution(A_current, B_current);
            
            // Paramètres LNS-v3 adaptés
            double neighborhood_ratio = 0.25;
            int max_iter_lns = 500;  // Moins d'itérations par alternance
            int stagnation_threshold = 10;
            
            // Temps alloué pour cette phase LNS (partager le temps restant)
            double lns_timeout = std::min(time_remaining_ms() * 0.4, 10000.0);  // Max 10s par phase LNS
            
            // LNS-v3 utilise déjà l'alternance Fix_A, Fix_B, Joint en interne
            LNSResult lns_result = bmf.solve_lns_v3(
                max_iter_lns,
                neighborhood_ratio,
                42 + total_alternations,
                stagnation_threshold,
                true,          // verbose (moins verbeux)
                lns_timeout,
                &g_interrupted
            );
            
            auto lns_end = std::chrono::high_resolution_clock::now();
            double lns_time = std::chrono::duration<double, std::milli>(lns_end - lns_start).count();
            v3_lns_time += lns_time;
            v3_iterations += lns_result.iterations;
            
            A_current = lns_result.A_solution;
            B_current = lns_result.B_solution;
            int errors_after_lns = lns_result.final_errors;
            
            std::cout << "  [LNS-v3] " << current_errors << " -> " << errors_after_lns 
                      << " erreurs (" << lns_result.iterations << " iters, " 
                      << std::fixed << std::setprecision(1) << lns_time << " ms)" << std::endl;
            
            current_errors = errors_after_lns;
            
            if (current_errors < best_errors) {
                best_errors = current_errors;
                A_best = A_current;
                B_best = B_current;
                stagnation_alternations = 0;
            } else {
                stagnation_alternations++;
            }
            
            if (current_errors == 0) {
                std::cout << "  Factorisation exacte trouvée par LNS-v3!" << std::endl;
                break;
            }
            
            // Vérifier la stagnation
            if (stagnation_alternations >= MAX_STAGNATION_ALTERNATIONS) {
                std::cout << "\n  Stagnation détectée après " << MAX_STAGNATION_ALTERNATIONS 
                          << " alternances sans amélioration." << std::endl;
                break;
            }
        }
        
        // Utiliser la meilleure solution trouvée
        v3_A_final = A_best;
        v3_B_final = B_best;
        v3_final_errors = best_errors;
        
        auto total_end = std::chrono::high_resolution_clock::now();
        v3_total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
        
        std::cout << "\n=== RÉSULTATS GREEDY + (WLS <-> LNS-v3)* ===" << std::endl;
        std::cout << "Alternances: " << total_alternations << std::endl;
        std::cout << "Erreurs: " << v3_init_errors << " -> " << v3_final_errors << std::endl;
        std::cout << "Itérations LNS-v3 totales: " << v3_iterations << std::endl;
        std::cout << "Temps WLS total: " << std::fixed << std::setprecision(1) << v3_ls_time << " ms" << std::endl;
        std::cout << "Temps LNS-v3 total: " << std::fixed << std::setprecision(1) << v3_lns_time << " ms" << std::endl;
        std::cout << "Temps total: " << std::fixed << std::setprecision(1) << v3_total_time << " ms" << std::endl;
        
        // Remplir le résultat
        result.k = k;
        result.init_errors = v3_init_errors;
        result.after_ls_errors = v3_init_errors;
        result.final_errors = v3_final_errors;
        result.iterations = v3_iterations;
        result.ls_time = v3_ls_time;
        result.lns_time = v3_lns_time;
        result.total_time = v3_total_time;
        result.success = true;
    }
    
    // ==================== COMPARAISON FINALE ====================
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== COMPARAISON FINALE ===" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "Fichier: " << CSVMatrixLoader::getBaseName(filename) << std::endl;
    std::cout << "Dimensions: " << m << "x" << n << ", k=" << k << std::endl;
    std::cout << "Alternances WLS<->LNS: " << total_alternations << std::endl;
    std::cout << std::endl;
    
    std::cout << std::setw(35) << "Méthode" 
              << std::setw(15) << "Init" 
              << std::setw(15) << "Final" 
              << std::setw(15) << "Temps (ms)" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    std::cout << std::setw(35) << "GREEDY + (WLS <-> LNS-v3)*" 
              << std::setw(15) << v3_init_errors 
              << std::setw(15) << v3_final_errors
              << std::setw(15) << std::fixed << std::setprecision(1) << v3_total_time << std::endl;
    
    std::cout << std::string(80, '-') << std::endl;
    
    std::cout << "\nRÉSULTAT: " << v3_final_errors << " erreurs finales" << std::endl;
    std::cout << "   Temps total: " << std::fixed << std::setprecision(1) << v3_total_time << " ms" << std::endl;
    
    // ==================== VÉRIFICATION PYTHON ====================
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== VÉRIFICATION PYTHON (verif.py) ===" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    int python_errors = verify_with_python(M, v3_A_final, v3_B_final, "GREEDY_WLS_LNS-v3_ALT", filename);
    
    if (python_errors >= 0) {
        if (python_errors == v3_final_errors) {
            std::cout << "Vérification OK: C++ (" << v3_final_errors << ") == Python (" << python_errors << ")" << std::endl;
        } else {
            std::cout << "Différence: C++ (" << v3_final_errors << ") != Python (" << python_errors << ")" << std::endl;
        }
    }
    
    return result;
}


// ==================== TEST MULTI-START ====================

/**
 * @brief Test la factorisation BMF avec multi-start (plusieurs seeds)
 * 
 * Lance N exécutions avec des seeds différents et retourne le meilleur résultat.
 * 
 * @param csv_file Chemin vers le fichier CSV
 * @param k_value Nombre de facteurs (0 = automatique)
 * @param num_starts Nombre d'exécutions avec seeds différents
 * @param timeout_seconds Timeout en secondes par exécution (0 = pas de timeout)
 * @return CSVTestResult contenant les résultats du meilleur run
 */
CSVTestResult test_csv_multistart(const std::string& csv_file, int k_value, int num_starts, double timeout_seconds = 0) {
    CSVTestResult best_result;
    best_result.filename = std::filesystem::path(csv_file).filename().string();
    best_result.k = k_value;
    best_result.success = false;
    best_result.final_errors = INT_MAX;
    
    // Installer le gestionnaire de signal pour Ctrl+C
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== MULTI-START: " << num_starts << " exécutions avec seeds différents ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::string filename = csv_file;
    
    // Si c'est juste un nom de fichier, ajouter le chemin data/
    if (filename.find('/') == std::string::npos) {
        filename = "../data/" + filename;
        if (filename.find(".csv") == std::string::npos) {
            filename += ".csv";
        }
    }
    
    std::cout << "Chargement de: " << filename << std::endl;
    
    Matrix M = CSVMatrixLoader::loadFromCSV(filename);
    DataGenerator::print_matrix_stats(M, "Matrice CSV");
    if (M.rows == 0 || M.cols == 0) {
        std::cerr << "Erreur: impossible de charger la matrice" << std::endl;
        return best_result;
    }
    
    int m = M.rows;
    int n = M.cols;
    
    // Déterminer k automatiquement si non spécifié
    int k = k_value;
    if (k <= 0) {
        k = std::max(5, std::min(50, static_cast<int>(std::log2(std::min(m, n)) * 2)));
        std::cout << "  k automatique: " << k << " (basé sur dimensions " << m << "x" << n << ")" << std::endl;
    } else {
        std::cout << "  k spécifié: " << k << std::endl;
    }
    
    if (timeout_seconds > 0) {
        std::cout << "  Timeout par run: " << std::fixed << std::setprecision(2) << timeout_seconds << " secondes" << std::endl;
    }
    
    // Stocker les résultats de chaque run
    std::vector<std::tuple<int, int, int, double>> run_results; // seed, wls_errors, final_errors, time
    
    int best_seed = 0;
    Matrix best_A(0, 0), best_B(0, 0);
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    for (int run = 0; run < num_starts && !is_interrupted(); run++) {
        int seed = 42 + run;  // Seeds: 42, 43, 44, ...
        
        std::cout << "\n" << std::string(50, '-') << std::endl;
        std::cout << ">>> RUN " << (run + 1) << "/" << num_starts << " (seed=" << seed << ")" << std::endl;
        std::cout << std::string(50, '-') << std::endl;
        
        auto run_start = std::chrono::high_resolution_clock::now();
        
        // Phase 1: GREEDY + WLS
        BMFLocalSearch ls(k, M, seed);
        ls.initialize_greedy();
        ls.compute_all_counts();
        int init_errors = ls.count_errors();
        
        int wls_errors = init_errors;
        if (init_errors > 0) {
            LocalSearchResult wls_result = ls.solve_weighted(
                10000000,   // max_iterations
                0.1,        // penalty_increment
                30,         // max_stagnation
                true,       // verbose (moins verbeux en multi-start)
                &g_interrupted  // stop_flag pour Ctrl+C
            );
            wls_errors = wls_result.final_errors;
            
            // IMPORTANT: Utiliser les MEILLEURES matrices trouvées par WLS
            ls.A = wls_result.A_solution;
            ls.B = wls_result.B_solution;
        }
        
        // Phase 2: LNS-v3
        int final_errors = wls_errors;
        Matrix A_final = ls.A;
        Matrix B_final = ls.B;
        
        if (wls_errors > 0) {
            BMF bmf(k, M);
            bmf.set_initial_solution(ls.A, ls.B);  // Maintenant ce sont les meilleures de WLS
            
            double timeout_ms = timeout_seconds * 1000.0;
            auto lns_result = bmf.solve_lns_v3(
                10000,     // max_iterations
                0.05,      // neighborhood_ratio
                seed,      // seed
                15,        // stagnation_threshold
                true,     // verbose
                timeout_ms,
                &g_interrupted
            );
            
            final_errors = lns_result.final_errors;
            A_final = lns_result.A_solution;
            B_final = lns_result.B_solution;
        }
        
        auto run_end = std::chrono::high_resolution_clock::now();
        double run_time = std::chrono::duration<double, std::milli>(run_end - run_start).count();
        
        run_results.push_back({seed, wls_errors, final_errors, run_time});
        
        std::cout << "  Résultat: " << init_errors << " -> " << wls_errors << " (WLS) -> " 
                  << final_errors << " erreurs en " << std::fixed << std::setprecision(1) 
                  << run_time << " ms" << std::endl;
        
        // Mettre à jour le meilleur
        if (final_errors < best_result.final_errors) {
            best_result.final_errors = final_errors;
            best_result.after_ls_errors = wls_errors;
            best_result.init_errors = init_errors;
            best_result.total_time = run_time;
            best_result.success = true;
            best_seed = seed;
            best_A = A_final;
            best_B = B_final;
            
            std::cout << "  NOUVEAU MEILLEUR! (seed=" << seed << ", " << final_errors << " erreurs)" << std::endl;
            
            // Si on atteint 0, on peut s'arrêter
            if (final_errors == 0) {
                std::cout << "  Solution parfaite trouvée, arrêt anticipé!" << std::endl;
                break;
            }
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // ==================== RÉSUMÉ MULTI-START ====================
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== RÉSUMÉ MULTI-START ===" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\n" << std::setw(10) << "Seed" << std::setw(15) << "Après WLS" 
              << std::setw(15) << "Final" << std::setw(15) << "Temps (ms)" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
    
    for (const auto& [seed, wls_err, final_err, time] : run_results) {
        std::string marker = (seed == best_seed) ? " " : "";
        std::cout << std::setw(10) << seed 
                  << std::setw(15) << wls_err 
                  << std::setw(15) << final_err 
                  << std::setw(15) << std::fixed << std::setprecision(1) << time 
                  << marker << std::endl;
    }
    
    std::cout << std::string(55, '-') << std::endl;
    std::cout << "\nMEILLEUR: seed=" << best_seed << " avec " << best_result.final_errors << " erreurs" << std::endl;
    std::cout << "   Temps total multi-start: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    
    // ==================== VÉRIFICATION PYTHON ====================
    if (best_result.success) {
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "=== VÉRIFICATION PYTHON (verif.py) ===" << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
        int python_errors = verify_with_python(M, best_A, best_B, "MULTI-START", filename);
        
        if (python_errors >= 0) {
            if (python_errors == best_result.final_errors) {
                std::cout << "Vérification OK: C++ (" << best_result.final_errors << ") == Python (" << python_errors << ")" << std::endl;
            } else {
                std::cout << "Différence: C++ (" << best_result.final_errors << ") != Python (" << python_errors << ")" << std::endl;
            }
        }
    }
    
    best_result.k = k;
    return best_result;
}


// ==================== TEST LS + WLS SEULEMENT ====================

/**
 * @brief Test la factorisation BMF avec GREEDY + WLS uniquement (sans LNS-v3)
 * 
 * @param csv_file Chemin vers le fichier CSV
 * @param k_value Nombre de facteurs (0 = automatique)
 * @return CSVTestResult contenant les résultats du test
 */
CSVTestResult test_csv_lswls(const std::string& csv_file, int k_value = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    // Installer le gestionnaire de signal pour Ctrl+C
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== GREEDY + WLS (sans LNS-v3) ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    std::cout << "  Pas de limite (Ctrl+C pour arreter)" << std::endl;
    
    std::string filename = csv_file;
    
    // Si c'est juste un nom de fichier, ajouter le chemin data/
    if (filename.find('/') == std::string::npos) {
        filename = "../data/" + filename;
        if (filename.find(".csv") == std::string::npos) {
            filename += ".csv";
        }
    }
    
    std::cout << "Chargement de: " << filename << std::endl;
    
    Matrix M = CSVMatrixLoader::loadFromCSV(filename);
    DataGenerator::print_matrix_stats(M, "Matrice CSV");
    if (M.rows == 0 || M.cols == 0) {
        std::cerr << "Erreur: impossible de charger la matrice" << std::endl;
        return result;
    }
    
    int m = M.rows;
    int n = M.cols;
    
    // Déterminer k automatiquement si non spécifié
    int k = k_value;
    if (k <= 0) {
        k = std::max(5, std::min(50, static_cast<int>(std::log2(std::min(m, n)) * 2)));
        std::cout << "  k automatique: " << k << " (basé sur dimensions " << m << "x" << n << ")" << std::endl;
    } else {
        std::cout << "  k spécifié: " << k << std::endl;
    }
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Phase 1: GREEDY
    std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
    
    BMFLocalSearch ls(k, M, 42);  // seed harmonisé
    ls.initialize_greedy();
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    
    std::cout << "  [GREEDY] Erreurs initiales: " << init_errors << std::endl;
    
    // Phase 2: WLS
    int wls_errors = init_errors;
    Matrix A_final = ls.A;
    Matrix B_final = ls.B;
    
    if (init_errors > 0) {
        std::cout << "\n[Phase 2] Local Search..." << std::endl;
        
        LocalSearchResult wls_result = ls.solve(
            5000,   // max_iterations
            true       // verbose
        );
        
        wls_errors = wls_result.final_errors;
        A_final = wls_result.A_solution;
        B_final = wls_result.B_solution;
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // Résultats
    std::cout << "\n=== RÉSULTATS LS + WLS ===" << std::endl;
    std::cout << "Erreurs: " << init_errors << " -> " << wls_errors << " (LS+WLS)" << std::endl;
    std::cout << "Temps total: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    
    result.k = k;
    result.init_errors = init_errors;
    result.after_ls_errors = wls_errors;
    result.final_errors = wls_errors;
    result.iterations = 0;  // Pas de LNS
    result.ls_time = total_time;
    result.lns_time = 0;
    result.total_time = total_time;
    result.success = true;
    
    // Vérification Python
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== VÉRIFICATION PYTHON ===" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    int python_errors = verify_with_python(M, A_final, B_final, "GREEDY_LS", filename);
    
    if (python_errors >= 0) {
        if (python_errors == wls_errors) {
            std::cout << " Vérification OK: C++ (" << wls_errors << ") == Python (" << python_errors << ")" << std::endl;
        } else {
            std::cout << "Différence: C++ (" << wls_errors << ") != Python (" << python_errors << ")" << std::endl;
        }
    }
    
    return result;
}


// ==================== TEST LS + LNS-v3 (sans WLS) ====================

/**
 * @brief Test la factorisation BMF avec GREEDY + LNS-v3 directement (sans WLS)
 * 
 * @param csv_file Chemin vers le fichier CSV
 * @param k_value Nombre de facteurs (0 = automatique)
 * @param timeout_seconds Timeout en secondes (0 = pas de timeout)
 * @return CSVTestResult contenant les résultats du test
 */
CSVTestResult test_csv_lslns(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    // Installer le gestionnaire de signal pour Ctrl+C
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== GREEDY + LNS-v3 (sans WLS) ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::string filename = csv_file;
    
    // Si c'est juste un nom de fichier, ajouter le chemin data/
    if (filename.find('/') == std::string::npos) {
        filename = "../data/" + filename;
        if (filename.find(".csv") == std::string::npos) {
            filename += ".csv";
        }
    }
    
    std::cout << "Chargement de: " << filename << std::endl;
    
    Matrix M = CSVMatrixLoader::loadFromCSV(filename);
    DataGenerator::print_matrix_stats(M, "Matrice CSV");
    if (M.rows == 0 || M.cols == 0) {
        std::cerr << "Erreur: impossible de charger la matrice" << std::endl;
        return result;
    }
    
    int m = M.rows;
    int n = M.cols;
    
    // Déterminer k automatiquement si non spécifié
    int k = k_value;
    if (k <= 0) {
        k = std::max(5, std::min(50, static_cast<int>(std::log2(std::min(m, n)) * 2)));
        std::cout << "  k automatique: " << k << " (basé sur dimensions " << m << "x" << n << ")" << std::endl;
    } else {
        std::cout << "  k spécifié: " << k << std::endl;
    }
    
    if (timeout_seconds > 0) {
        std::cout << "  Timeout LNS-v3: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "  Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Phase 1: GREEDY seulement
    std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
    auto ls_start = std::chrono::high_resolution_clock::now();
    
    BMFLocalSearch ls(k, M, 43);  // seed harmonisé
    ls.initialize_greedy();
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    
    auto ls_end = std::chrono::high_resolution_clock::now();
    double ls_time = std::chrono::duration<double, std::milli>(ls_end - ls_start).count();
    
    std::cout << "  [GREEDY] Erreurs initiales: " << init_errors << " (" << std::fixed << std::setprecision(1) << ls_time << " ms)" << std::endl;
    
    // Phase 2: LNS-v3 directement (sans passer par WLS)
    int final_errors = init_errors;
    int iterations = 0;
    double lns_time = 0;
    Matrix A_final = ls.A;
    Matrix B_final = ls.B;
    
    if (init_errors > 0) {
        std::cout << "\n[Phase 2] LNS-v3 (Alternance A/B + Anti-stagnation)..." << std::endl;
        
        auto lns_start = std::chrono::high_resolution_clock::now();
        
        BMF bmf(k, M);
        bmf.set_initial_solution(ls.A, ls.B);  // Solution GREEDY directement
        
        // Paramètres LNS-v3
        double neighborhood_ratio = 0.05;  // 5% des dimensions
        int max_iter_v3 = 10000;
        int stagnation_threshold = 15;
        double timeout_ms = timeout_seconds * 1000.0;
        
        auto lns_result = bmf.solve_lns_v3(
            max_iter_v3,
            neighborhood_ratio,
            43,               // seed
            stagnation_threshold,
            true,             // verbose
            timeout_ms,
            &g_interrupted
        );
        
        auto lns_end = std::chrono::high_resolution_clock::now();
        lns_time = std::chrono::duration<double, std::milli>(lns_end - lns_start).count();
        final_errors = lns_result.final_errors;
        iterations = lns_result.iterations;
        A_final = lns_result.A_solution;
        B_final = lns_result.B_solution;
        
        // Afficher la raison de l'arrêt
        std::string stop_label = (lns_result.stop_reason == "zero_errors") ? "[OK]" :
                                 (lns_result.stop_reason == "timeout") ? "[TIMEOUT]" :
                                 (lns_result.stop_reason == "interrupted") ? "[INTERROMPU]" : "[STAGNATION]";
        
        std::cout << "  LNS-v3: " << init_errors << " -> " << final_errors 
                  << " erreurs en " << iterations << " iters ("
                  << std::fixed << std::setprecision(1) << lns_time << " ms)" << std::endl;
        std::cout << "  Arret: " << stop_label << " " << lns_result.stop_reason << std::endl;
    } else {
        std::cout << "\n[Phase 2] Pas besoin de LNS-v3, déjà 0 erreurs!" << std::endl;
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // Résultats
    std::cout << "\n=== RÉSULTATS LS + LNS-v3 ===" << std::endl;
    std::cout << "Erreurs: " << init_errors << " (LS) -> " << final_errors << " (LNS-v3)" << std::endl;
    std::cout << "Itérations LNS-v3: " << iterations << std::endl;
    std::cout << "Temps LS: " << std::fixed << std::setprecision(1) << ls_time << " ms" << std::endl;
    std::cout << "Temps LNS-v3: " << std::fixed << std::setprecision(1) << lns_time << " ms" << std::endl;
    std::cout << "Temps total: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    
    result.k = k;
    result.init_errors = init_errors;
    result.after_ls_errors = init_errors;  // Pas de WLS
    result.final_errors = final_errors;
    result.iterations = iterations;
    result.ls_time = ls_time;
    result.lns_time = lns_time;
    result.total_time = total_time;
    result.success = true;
    
    // Vérification Python
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== VÉRIFICATION PYTHON ===" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    int python_errors = verify_with_python(M, A_final, B_final, "GREEDY_LNS-v3", filename);
    
    if (python_errors >= 0) {
        if (python_errors == final_errors) {
            std::cout << "Vérification OK: C++ (" << final_errors << ") == Python (" << python_errors << ")" << std::endl;
        } else {
            std::cout << "Différence: C++ (" << final_errors << ") != Python (" << python_errors << ")" << std::endl;
        }
    }
    
    return result;
}


// ==================== OPTIONS MODULAIRES 1-6 ====================

/**
 * Helper: Charge un fichier CSV et retourne la matrice M avec les infos
 */
std::tuple<Matrix, int, std::string> load_csv_matrix(const std::string& csv_file, int k_value) {
    std::string filename = csv_file;
    
    // Si c'est juste un nom de fichier, ajouter le chemin data/
    if (filename.find('/') == std::string::npos) {
        filename = "../data/" + filename;
        if (filename.find(".csv") == std::string::npos) {
            filename += ".csv";
        }
    }
    
    std::cout << "Chargement de: " << filename << std::endl;
    
    Matrix M = CSVMatrixLoader::loadFromCSV(filename);
    DataGenerator::print_matrix_stats(M, "Matrice CSV");
    
    if (M.rows == 0 || M.cols == 0) {
        std::cerr << "Erreur: impossible de charger la matrice" << std::endl;
        return {Matrix(0, 0), 0, filename};
    }
    
    // Déterminer k automatiquement si non spécifié
    int k = k_value;
    if (k <= 0) {
        k = std::max(5, std::min(50, static_cast<int>(std::log2(std::min(M.rows, M.cols)) * 2)));
        std::cout << "  k automatique: " << k << " (basé sur dimensions " << M.rows << "x" << M.cols << ")" << std::endl;
    } else {
        std::cout << "  k spécifié: " << k << std::endl;
    }
    
    return {M, k, filename};
}

/**
 * Helper: Sélectionne des lignes/colonnes pour LNS avec stratégie 70% top-k + 30% random
 * @param M La matrice cible
 * @param A_current Solution actuelle A
 * @param B_current Solution actuelle B
 * @param ratio Proportion de lignes/colonnes à sélectionner
 * @param rng Générateur aléatoire
 * @return Paire (lignes sélectionnées, colonnes sélectionnées)
 */
std::pair<std::vector<int>, std::vector<int>> select_neighborhood_topk_random(
    const Matrix& M, const Matrix& A_current, const Matrix& B_current,
    int count, std::mt19937& rng
) {
    int m = M.rows;
    int n = M.cols;
    int num_rows = std::min(count, m);
    int num_cols = std::min(count, n);
    
    // Calculer les erreurs par ligne et par colonne
    Matrix computed = A_current.multiply(B_current);
    std::vector<std::pair<int, int>> row_errors(m);  // (erreurs, index)
    std::vector<std::pair<int, int>> col_errors(n);  // (erreurs, index)
    
    for (int i = 0; i < m; i++) {
        int errors = 0;
        for (int j = 0; j < n; j++) {
            if (M(i, j) != -1 && computed(i, j) != M(i, j)) {
                errors++;
            }
        }
        row_errors[i] = {errors, i};
    }
    
    for (int j = 0; j < n; j++) {
        int errors = 0;
        for (int i = 0; i < m; i++) {
            if (M(i, j) != -1 && computed(i, j) != M(i, j)) {
                errors++;
            }
        }
        col_errors[j] = {errors, j};
    }
    
    // Trier par nombre d'erreurs décroissant
    std::sort(row_errors.begin(), row_errors.end(), std::greater<>());
    std::sort(col_errors.begin(), col_errors.end(), std::greater<>());
    
    // 70% top-k, 30% random
    int topk_rows = std::max(1, static_cast<int>(num_rows * 0.7));
    int random_rows = num_rows - topk_rows;
    int topk_cols = std::max(1, static_cast<int>(num_cols * 0.7));
    int random_cols = num_cols - topk_cols;
    
    std::set<int> selected_rows_set;
    std::set<int> selected_cols_set;
    
    // Ajouter les top-k lignes (celles avec le plus d'erreurs)
    for (int i = 0; i < topk_rows && i < m; i++) {
        selected_rows_set.insert(row_errors[i].second);
    }
    
    // Ajouter les top-k colonnes
    for (int j = 0; j < topk_cols && j < n; j++) {
        selected_cols_set.insert(col_errors[j].second);
    }
    
    // Ajouter des lignes aléatoires (parmi celles non encore sélectionnées)
    std::vector<int> remaining_rows;
    for (int i = topk_rows; i < m; i++) {
        remaining_rows.push_back(row_errors[i].second);
    }
    std::shuffle(remaining_rows.begin(), remaining_rows.end(), rng);
    for (int i = 0; i < random_rows && i < static_cast<int>(remaining_rows.size()); i++) {
        selected_rows_set.insert(remaining_rows[i]);
    }
    
    // Ajouter des colonnes aléatoires
    std::vector<int> remaining_cols;
    for (int j = topk_cols; j < n; j++) {
        remaining_cols.push_back(col_errors[j].second);
    }
    std::shuffle(remaining_cols.begin(), remaining_cols.end(), rng);
    for (int j = 0; j < random_cols && j < static_cast<int>(remaining_cols.size()); j++) {
        selected_cols_set.insert(remaining_cols[j]);
    }
    
    std::vector<int> selected_rows(selected_rows_set.begin(), selected_rows_set.end());
    std::vector<int> selected_cols(selected_cols_set.begin(), selected_cols_set.end());
    
    return {selected_rows, selected_cols};
}

/**
 * Helper: Sélectionne des COLONNES avec stratégie 70% top-k + 30% random (pour Fix_A)
 * Les colonnes avec le plus d'erreurs sont prioritaires.
 */
std::vector<int> select_cols_topk_random(
    const Matrix& M, const Matrix& A_current, const Matrix& B_current,
    int count, std::mt19937& rng
) {
    int m = M.rows;
    int n = M.cols;
    int num_cols = std::min(count, n);
    
    // Calculer les erreurs par colonne
    Matrix computed = A_current.multiply(B_current);
    std::vector<std::pair<int, int>> col_errors(n);  // (erreurs, index)
    
    for (int j = 0; j < n; j++) {
        int errors = 0;
        for (int i = 0; i < m; i++) {
            if (M(i, j) != -1 && computed(i, j) != M(i, j)) {
                errors++;
            }
        }
        col_errors[j] = {errors, j};
    }
    
    // Trier par nombre d'erreurs décroissant
    std::sort(col_errors.begin(), col_errors.end(), std::greater<>());
    
    // 70% top-k, 30% random
    int topk_cols = std::max(1, static_cast<int>(num_cols * 0.7));
    int random_cols = num_cols - topk_cols;
    
    std::set<int> selected_set;
    
    // Ajouter les top-k colonnes (celles avec le plus d'erreurs)
    for (int j = 0; j < topk_cols && j < n; j++) {
        selected_set.insert(col_errors[j].second);
    }
    
    // Ajouter des colonnes aléatoires (parmi les restantes)
    std::vector<int> remaining;
    for (int j = topk_cols; j < n; j++) {
        remaining.push_back(col_errors[j].second);
    }
    std::shuffle(remaining.begin(), remaining.end(), rng);
    for (int j = 0; j < random_cols && j < static_cast<int>(remaining.size()); j++) {
        selected_set.insert(remaining[j]);
    }
    
    return std::vector<int>(selected_set.begin(), selected_set.end());
}

/**
 * Helper: Sélectionne des LIGNES avec stratégie 70% top-k + 30% random (pour Fix_B)
 * Les lignes avec le plus d'erreurs sont prioritaires.
 */
std::vector<int> select_rows_topk_random(
    const Matrix& M, const Matrix& A_current, const Matrix& B_current,
    int count, std::mt19937& rng
) {
    int m = M.rows;
    int n = M.cols;
    int num_rows = std::min(count, m);
    
    // Calculer les erreurs par ligne
    Matrix computed = A_current.multiply(B_current);
    std::vector<std::pair<int, int>> row_errors(m);  // (erreurs, index)
    
    for (int i = 0; i < m; i++) {
        int errors = 0;
        for (int j = 0; j < n; j++) {
            if (M(i, j) != -1 && computed(i, j) != M(i, j)) {
                errors++;
            }
        }
        row_errors[i] = {errors, i};
    }
    
    // Trier par nombre d'erreurs décroissant
    std::sort(row_errors.begin(), row_errors.end(), std::greater<>());
    
    // 70% top-k, 30% random
    int topk_rows = std::max(1, static_cast<int>(num_rows * 0.7));
    int random_rows = num_rows - topk_rows;
    
    std::set<int> selected_set;
    
    // Ajouter les top-k lignes (celles avec le plus d'erreurs)
    for (int i = 0; i < topk_rows && i < m; i++) {
        selected_set.insert(row_errors[i].second);
    }
    
    // Ajouter des lignes aléatoires (parmi les restantes)
    std::vector<int> remaining;
    for (int i = topk_rows; i < m; i++) {
        remaining.push_back(row_errors[i].second);
    }
    std::shuffle(remaining.begin(), remaining.end(), rng);
    for (int i = 0; i < random_rows && i < static_cast<int>(remaining.size()); i++) {
        selected_set.insert(remaining[i]);
    }
    
    return std::vector<int>(selected_set.begin(), selected_set.end());
}

/**
 * Helper: Compte les erreurs entre M et A×B
 */
int count_reconstruction_errors(const Matrix& M, const Matrix& A, const Matrix& B) {
    Matrix computed = A.multiply(B);
    int errors = 0;
    for (int i = 0; i < M.rows; i++) {
        for (int j = 0; j < M.cols; j++) {
            if (M(i, j) != -1 && computed(i, j) != M(i, j)) {
                errors++;
            }
        }
    }
    return errors;
}

/**
 * @brief Option 1: GREEDY + Local Search Basic (jusqu'au minimum local)
 */
CSVTestResult test_opt1_greedy_basic(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== OPTION 1: GREEDY + Local Search BASIC ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    auto [M, k, filename] = load_csv_matrix(csv_file, k_value);
    if (M.rows == 0) return result;
    
    if (timeout_seconds > 0) {
        std::cout << "  Timeout: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "  Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Phase 1: GREEDY
    std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
    BMFLocalSearch ls(k, M, 42);
    ls.initialize_greedy();
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    std::cout << "  [GREEDY] Erreurs: " << init_errors << std::endl;
    
    // Phase 2: Local Search BASIC (jusqu'au minimum local)
    int final_errors = init_errors;
    Matrix A_final = ls.A;
    Matrix B_final = ls.B;
    
    if (init_errors > 0 && !is_interrupted()) {
        std::cout << "\n[Phase 2] Local Search BASIC (jusqu'au minimum local)..." << std::endl;
        
        LocalSearchResult basic_result = ls.solve_basic(0, true);  // 0 = pas de limite
        
        final_errors = basic_result.final_errors;
        A_final = basic_result.A_solution;
        B_final = basic_result.B_solution;
        
        std::cout << "  [BASIC] " << init_errors << " -> " << final_errors 
                  << " erreurs en " << basic_result.iterations << " iters ("
                  << std::fixed << std::setprecision(1) << basic_result.total_time << " ms)" << std::endl;
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // Résultats
    std::cout << "\n=== RÉSULTATS OPTION 1 ===" << std::endl;
    std::cout << "Erreurs: " << init_errors << " -> " << final_errors << std::endl;
    std::cout << "Temps total: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    
    result.k = k;
    result.init_errors = init_errors;
    result.final_errors = final_errors;
    result.total_time = total_time;
    result.success = true;
    
    // Vérification Python
    int python_errors = verify_with_python(M, A_final, B_final, "OPT1_GREEDY_BASIC", filename);
    if (python_errors >= 0 && python_errors == final_errors) {
        std::cout << "Vérification OK" << std::endl;
    }
    
    return result;
}

/**
 * @brief Option 2: GREEDY + Fix_A/Fix_B (lns_step_fix_A et lns_step_fix_B alternés)
 */
CSVTestResult test_opt2_greedy_fix_ab(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== OPTION 2: GREEDY + Fix_A/Fix_B (lns_step_fix_A/fix_B) ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    auto [M, k, filename] = load_csv_matrix(csv_file, k_value);
    if (M.rows == 0) return result;
    
    auto total_start = std::chrono::high_resolution_clock::now();
    double timeout_ms = timeout_seconds * 1000.0;
    
    auto time_remaining_ms = [&]() -> double {
        if (timeout_ms <= 0) return std::numeric_limits<double>::max();
        auto now = std::chrono::high_resolution_clock::now();
        return timeout_ms - std::chrono::duration<double, std::milli>(now - total_start).count();
    };
    
    if (timeout_seconds > 0) {
        std::cout << "  Timeout: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "  Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    // Phase 1: GREEDY
    std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
    BMFLocalSearch ls(k, M, 42);
    ls.initialize_greedy();
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    std::cout << "  [GREEDY] Erreurs: " << init_errors << std::endl;
    
    // Phase 2: Fix_A/Fix_B alternés
    int best_errors = init_errors;
    Matrix A_best = ls.A;
    Matrix B_best = ls.B;
    Matrix A_current = ls.A;
    Matrix B_current = ls.B;
    int total_iterations = 0;
    std::mt19937 rng(42);
    int neighborhood_size = 10;  // taille du voisinage (10, 20, 30, ..., 100)
    int max_neighborhood = std::min(100, std::min(M.rows, M.cols));  // plafonné à 100 ou dimension min
    std::cout << "  Voisinage: " << neighborhood_size << " -> " << max_neighborhood 
              << " (par pas de 10, matrice " << M.rows << "x" << M.cols << ")" << std::endl;
    
    BMF bmf(k, M);
    
    while (best_errors > 0 && !is_interrupted() && time_remaining_ms() > 1000) {
        total_iterations++;
        bool improved = false;
        
        // Sélectionner des colonnes (70% top-k erreurs + 30% random) pour Fix_A
        std::vector<int> selected_cols = select_cols_topk_random(M, A_current, B_current, neighborhood_size, rng);
        
        // Fix_A PROGRESSIF: optimiser B colonne par colonne (A fixé)
        std::cout << "[Fix_A] Iter " << total_iterations << ": " << selected_cols.size() << " cols (voisinage=" << neighborhood_size << ")" << std::endl;
        for (int col_j : selected_cols) {
            std::vector<int> single_col = {col_j};
            int new_cost = bmf.lns_step_fix_A(A_current, B_current, single_col);
            if (new_cost >= 0 && new_cost < best_errors) {
                best_errors = new_cost;
                A_best = A_current;
                B_best = B_current;
                improved = true;
                std::cout << "  \U0001F3AF Fix_A col " << col_j << ": " << best_errors << " erreurs" << std::endl;
                if (best_errors == 0) break;
            }
        }
        
        if (best_errors == 0 || is_interrupted() || time_remaining_ms() < 1000) break;
        
        // Sélectionner des lignes (70% top-k erreurs + 30% random) pour Fix_B
        std::vector<int> selected_rows = select_rows_topk_random(M, A_current, B_current, neighborhood_size, rng);
        
        // Fix_B PROGRESSIF: optimiser A ligne par ligne (B fixé)
        std::cout << "[Fix_B] Iter " << total_iterations << ": " << selected_rows.size() << " rows (voisinage=" << neighborhood_size << ")" << std::endl;
        for (int row_i : selected_rows) {
            std::vector<int> single_row = {row_i};
            int new_cost = bmf.lns_step_fix_B(A_current, B_current, single_row);
            if (new_cost >= 0 && new_cost < best_errors) {
                best_errors = new_cost;
                A_best = A_current;
                B_best = B_current;
                improved = true;
                std::cout << "  \U0001F3AF Fix_B row " << row_i << ": " << best_errors << " erreurs" << std::endl;
                if (best_errors == 0) break;
            }
        }
        
        // Détection de minimum local
        if (!improved) {
            if (neighborhood_size >= max_neighborhood) {
                std::cout << "  Minimum local atteint à iter " << total_iterations 
                          << " avec " << best_errors << " erreurs (voisinage=" << max_neighborhood << " = max)" << std::endl;
                break;
            }
            int old_size = neighborhood_size;
            neighborhood_size = std::min(max_neighborhood, neighborhood_size + 10);
            std::cout << "  Pas d'amélioration, voisinage élargi: " 
                      << old_size << " -> " << neighborhood_size << std::endl;
        } else {
            neighborhood_size = 10;  // reset après amélioration
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // Résultats
    std::cout << "\n=== RÉSULTATS OPTION 2 ===" << std::endl;
    std::cout << "Erreurs: " << init_errors << " -> " << best_errors << std::endl;
    std::cout << "Itérations: " << total_iterations << std::endl;
    std::cout << "Temps total: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    
    result.k = k;
    result.init_errors = init_errors;
    result.final_errors = best_errors;
    result.iterations = total_iterations;
    result.total_time = total_time;
    result.success = true;
    
    // Vérification Python
    int python_errors = verify_with_python(M, A_best, B_best, "OPT2_GREEDY_FIX_AB", filename);
    if (python_errors >= 0 && python_errors == best_errors) {
        std::cout << "Vérification OK" << std::endl;
    }
    
    return result;
}

/**
 * @brief Option 3: GREEDY + (LNS)* en boucle (lns_step_partial uniquement)
 */
CSVTestResult test_opt3_greedy_lns_loop(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== OPTION 3: GREEDY + (LNS)* en boucle (lns_step_partial) ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    auto [M, k, filename] = load_csv_matrix(csv_file, k_value);
    if (M.rows == 0) return result;
    
    auto total_start = std::chrono::high_resolution_clock::now();
    double timeout_ms = timeout_seconds * 1000.0;
    
    auto time_remaining_ms = [&]() -> double {
        if (timeout_ms <= 0) return std::numeric_limits<double>::max();
        auto now = std::chrono::high_resolution_clock::now();
        return timeout_ms - std::chrono::duration<double, std::milli>(now - total_start).count();
    };
    
    if (timeout_seconds > 0) {
        std::cout << "  Timeout: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "  Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    // Phase 1: GREEDY
    std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
    BMFLocalSearch ls(k, M, 42);
    ls.initialize_greedy();
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    std::cout << "  [GREEDY] Erreurs: " << init_errors << std::endl;
    
    // Boucle LNS (lns_step_partial)
    int best_errors = init_errors;
    Matrix A_best = ls.A;
    Matrix B_best = ls.B;
    Matrix A_current = ls.A;
    Matrix B_current = ls.B;
    int total_iterations = 0;
    int stagnation = 0;
    std::mt19937 rng(42);
    int neighborhood_size = 10;  // taille du voisinage (10, 20, ..., 100)
    int max_neighborhood = std::min(100, std::min(M.rows, M.cols));
    
    BMF bmf(k, M);
    
    while (best_errors > 0 && !is_interrupted() && time_remaining_ms() > 1000) {
        total_iterations++;
        
        // Sélectionner un voisinage 70% top-k + 30% random
        auto [selected_rows, selected_cols] = select_neighborhood_topk_random(M, A_current, B_current, neighborhood_size, rng);
        
        std::cout << "[LNS] Iter " << total_iterations << ": partial avec " 
                  << selected_rows.size() << " lignes, " << selected_cols.size() << " colonnes..." << std::endl;
        
        int new_cost = bmf.lns_step_partial(A_current, B_current, selected_rows, selected_cols);
        
        if (new_cost >= 0 && new_cost < best_errors) {
            best_errors = new_cost;
            A_best = A_current;
            B_best = B_current;
            stagnation = 0;
            std::cout << "  Nouveau meilleur: " << best_errors << " erreurs" << std::endl;
        } else {
            stagnation++;
            // Augmenter le voisinage si stagnation
            if (stagnation >= 5) {
                neighborhood_size = std::min(max_neighborhood, neighborhood_size + 10);
                stagnation = 0;
            }
        }
        
        // Utiliser la meilleure solution
        A_current = A_best;
        B_current = B_best;
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // Résultats
    std::cout << "\n=== RÉSULTATS OPTION 3 ===" << std::endl;
    std::cout << "Erreurs: " << init_errors << " -> " << best_errors << std::endl;
    std::cout << "Itérations LNS: " << total_iterations << std::endl;
    std::cout << "Temps total: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    
    result.k = k;
    result.init_errors = init_errors;
    result.final_errors = best_errors;
    result.iterations = total_iterations;
    result.total_time = total_time;
    result.success = true;
    
    // Vérification Python
    int python_errors = verify_with_python(M, A_best, B_best, "OPT3_GREEDY_LNS_LOOP", filename);
    if (python_errors >= 0 && python_errors == best_errors) {
        std::cout << "Vérification OK" << std::endl;
    }
    
    return result;
}

/**
 * @brief Option 4: GREEDY + (Basic + LNS)* en boucle
 * Basic = solve_basic, LNS = lns_step_partial
 */
CSVTestResult test_opt4_greedy_basic_lns_loop(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== OPTION 4: GREEDY + (Basic + LNS)* en boucle ===" << std::endl;
    std::cout << "=== Basic = solve_basic, LNS = lns_step_partial ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    auto [M, k, filename] = load_csv_matrix(csv_file, k_value);
    if (M.rows == 0) return result;
    
    auto total_start = std::chrono::high_resolution_clock::now();
    double timeout_ms = timeout_seconds * 1000.0;
    
    auto time_remaining_ms = [&]() -> double {
        if (timeout_ms <= 0) return std::numeric_limits<double>::max();
        auto now = std::chrono::high_resolution_clock::now();
        return timeout_ms - std::chrono::duration<double, std::milli>(now - total_start).count();
    };
    
    if (timeout_seconds > 0) {
        std::cout << "  Timeout: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "  Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    // Phase 1: GREEDY
    std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
    BMFLocalSearch ls(k, M, 42);
    ls.initialize_greedy();
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    std::cout << "  [GREEDY] Erreurs: " << init_errors << std::endl;
    
    // Boucle (Basic + LNS)*
    int best_errors = init_errors;
    Matrix A_best = ls.A;
    Matrix B_best = ls.B;
    Matrix A_current = ls.A;
    Matrix B_current = ls.B;
    int total_loops = 0;
    int total_lns_iterations = 0;
    std::mt19937 rng(42);
    int neighborhood_size = 10;
    int max_neighborhood = std::min(100, std::min(M.rows, M.cols));
    
    BMF bmf(k, M);
    
    while (best_errors > 0 && !is_interrupted() && time_remaining_ms() > 1000) {
        total_loops++;
        
        std::cout << "\n=== Boucle " << total_loops << " (erreurs: " << best_errors << ") ===" << std::endl;
        
        // Phase A: Local Search BASIC (solve_basic)
        std::cout << "[BASIC] Local Search jusqu'au minimum local..." << std::endl;
        BMFLocalSearch ls_basic(k, M, 42 + total_loops);
        ls_basic.A = A_current;
        ls_basic.B = B_current;
        ls_basic.compute_all_counts();
        
        LocalSearchResult basic_result = ls_basic.solve_basic(0, true);
        
        if (basic_result.final_errors < best_errors) {
            best_errors = basic_result.final_errors;
            A_best = basic_result.A_solution;
            B_best = basic_result.B_solution;
            std::cout << "  Nouveau meilleur après BASIC: " << best_errors << " erreurs" << std::endl;
        }
        
        A_current = basic_result.A_solution;
        B_current = basic_result.B_solution;
        
        if (best_errors == 0) break;
        if (is_interrupted() || time_remaining_ms() < 1000) break;
        
        // Phase B: LNS (lns_step_partial)
        std::cout << "[LNS] lns_step_partial..." << std::endl;
        
        auto [selected_rows, selected_cols] = select_neighborhood_topk_random(M, A_current, B_current, neighborhood_size, rng);
        
        int new_cost = bmf.lns_step_partial(A_current, B_current, selected_rows, selected_cols);
        total_lns_iterations++;
        
        if (new_cost >= 0 && new_cost < best_errors) {
            best_errors = new_cost;
            A_best = A_current;
            B_best = B_current;
            std::cout << "  Nouveau meilleur après LNS: " << best_errors << " erreurs" << std::endl;
        } else {
            // Augmenter le voisinage si pas d'amélioration
            neighborhood_size = std::min(max_neighborhood, neighborhood_size + 10);
        }
        
        A_current = A_best;
        B_current = B_best;
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // Résultats
    std::cout << "\n=== RÉSULTATS OPTION 4 ===" << std::endl;
    std::cout << "Erreurs: " << init_errors << " -> " << best_errors << std::endl;
    std::cout << "Boucles: " << total_loops << std::endl;
    std::cout << "Itérations LNS: " << total_lns_iterations << std::endl;
    std::cout << "Temps total: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    
    result.k = k;
    result.init_errors = init_errors;
    result.final_errors = best_errors;
    result.iterations = total_lns_iterations;
    result.total_time = total_time;
    result.success = true;
    
    int python_errors = verify_with_python(M, A_best, B_best, "OPT4_GREEDY_BASIC_LNS", filename);
    if (python_errors >= 0 && python_errors == best_errors) {
        std::cout << "Vérification OK" << std::endl;
    }
    
    return result;
}

/**
 * @brief Option 5: GREEDY + (Fix_A/Fix_B + LNS)* en boucle
 * Fix_A/Fix_B = lns_step_fix_A/fix_B, LNS = lns_step_partial
 */
CSVTestResult test_opt5_greedy_fix_lns_loop(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== OPTION 5: GREEDY + (Fix_A/Fix_B + LNS)* en boucle ===" << std::endl;
    std::cout << "=== Fix_A/B = lns_step_fix_A/B, LNS = lns_step_partial ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    auto [M, k, filename] = load_csv_matrix(csv_file, k_value);
    if (M.rows == 0) return result;
    
    auto total_start = std::chrono::high_resolution_clock::now();
    double timeout_ms = timeout_seconds * 1000.0;
    
    auto time_remaining_ms = [&]() -> double {
        if (timeout_ms <= 0) return std::numeric_limits<double>::max();
        auto now = std::chrono::high_resolution_clock::now();
        return timeout_ms - std::chrono::duration<double, std::milli>(now - total_start).count();
    };
    
    if (timeout_seconds > 0) {
        std::cout << "  Timeout: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "  Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    // Phase 1: GREEDY
    std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
    BMFLocalSearch ls(k, M, 42);
    ls.initialize_greedy();
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    std::cout << "  [GREEDY] Erreurs: " << init_errors << std::endl;
    
    // Boucle (Fix_A/Fix_B + LNS)*
    int best_errors = init_errors;
    Matrix A_best = ls.A;
    Matrix B_best = ls.B;
    Matrix A_current = ls.A;
    Matrix B_current = ls.B;
    int total_loops = 0;
    int total_iterations = 0;
    std::mt19937 rng(42);
    int neighborhood_size = 10;
    int max_neighborhood = std::min(100, std::min(M.rows, M.cols));
    std::cout << "  Voisinage: " << neighborhood_size << " -> " << max_neighborhood 
              << " (par pas de 10, matrice " << M.rows << "x" << M.cols << ")" << std::endl;
    
    BMF bmf(k, M);
    
    while (best_errors > 0 && !is_interrupted() && time_remaining_ms() > 1000) {
        total_loops++;
        bool loop_improved = false;
        
        std::cout << "\n=== Boucle " << total_loops << " (erreurs: " << best_errors << ", voisinage=" << neighborhood_size << ") ===" << std::endl;
        
        // Phase A: Fix_A PROGRESSIF (colonne par colonne) - 70% top-k + 30% random
        std::vector<int> selected_cols = select_cols_topk_random(M, A_current, B_current, neighborhood_size, rng);
        std::cout << "[Fix_A] " << selected_cols.size() << " cols..." << std::endl;
        
        for (int col_j : selected_cols) {
            std::vector<int> single_col = {col_j};
            int new_cost = bmf.lns_step_fix_A(A_current, B_current, single_col);
            total_iterations++;
            if (new_cost >= 0 && new_cost < best_errors) {
                best_errors = new_cost;
                A_best = A_current;
                B_best = B_current;
                loop_improved = true;
                std::cout << "  \U0001F3AF Fix_A col " << col_j << ": " << best_errors << " erreurs" << std::endl;
                if (best_errors == 0) break;
            }
        }
        
        if (best_errors == 0) break;
        if (is_interrupted() || time_remaining_ms() < 1000) break;
        
        // Phase B: Fix_B PROGRESSIF (ligne par ligne) - 70% top-k + 30% random
        std::vector<int> selected_rows = select_rows_topk_random(M, A_current, B_current, neighborhood_size, rng);
        std::cout << "[Fix_B] " << selected_rows.size() << " rows..." << std::endl;
        
        for (int row_i : selected_rows) {
            std::vector<int> single_row = {row_i};
            int new_cost = bmf.lns_step_fix_B(A_current, B_current, single_row);
            total_iterations++;
            if (new_cost >= 0 && new_cost < best_errors) {
                best_errors = new_cost;
                A_best = A_current;
                B_best = B_current;
                loop_improved = true;
                std::cout << "  \U0001F3AF Fix_B row " << row_i << ": " << best_errors << " erreurs" << std::endl;
                if (best_errors == 0) break;
            }
        }
        
        if (best_errors == 0) break;
        if (is_interrupted() || time_remaining_ms() < 1000) break;
        
        // Phase C: LNS (lns_step_partial)
        std::cout << "[LNS] lns_step_partial (" << neighborhood_size << "x" << neighborhood_size << ")..." << std::endl;
        auto [rows_free, cols_free] = select_neighborhood_topk_random(M, A_current, B_current, neighborhood_size, rng);
        
        int new_cost = bmf.lns_step_partial(A_current, B_current, rows_free, cols_free);
        total_iterations++;
        
        if (new_cost >= 0 && new_cost < best_errors) {
            best_errors = new_cost;
            A_best = A_current;
            B_best = B_current;
            loop_improved = true;
            std::cout << "  LNS: " << best_errors << " erreurs" << std::endl;
        }
        
        A_current = A_best;
        B_current = B_best;
        
        // Détection de minimum local
        if (!loop_improved) {
            if (neighborhood_size >= max_neighborhood) {
                std::cout << "  Minimum local atteint à boucle " << total_loops 
                          << " avec " << best_errors << " erreurs (voisinage=" << max_neighborhood << " = max)" << std::endl;
                break;
            }
            int old_size = neighborhood_size;
            neighborhood_size = std::min(max_neighborhood, neighborhood_size + 10);
            std::cout << "  Pas d'amélioration, voisinage élargi: " 
                      << old_size << " -> " << neighborhood_size << std::endl;
        } else {
            neighborhood_size = 10;
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // Résultats
    std::cout << "\n=== RÉSULTATS OPTION 5 ===" << std::endl;
    std::cout << "Erreurs: " << init_errors << " -> " << best_errors << std::endl;
    std::cout << "Boucles: " << total_loops << std::endl;
    std::cout << "Itérations totales: " << total_iterations << std::endl;
    std::cout << "Temps total: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    
    result.k = k;
    result.init_errors = init_errors;
    result.final_errors = best_errors;
    result.iterations = total_iterations;
    result.total_time = total_time;
    result.success = true;
    
    int python_errors = verify_with_python(M, A_best, B_best, "OPT5_GREEDY_FIX_LNS", filename);
    if (python_errors >= 0 && python_errors == best_errors) {
        std::cout << "Vérification OK" << std::endl;
    }
    
    return result;
}

/**
 * @brief Option 6: GREEDY + (Basic + Fix_A/Fix_B + LNS)* en boucle
 * Basic = solve_basic, Fix_A/Fix_B = lns_step_fix_A/fix_B, LNS = lns_step_partial
 */
CSVTestResult test_opt6_greedy_full_loop(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== OPTION 6: GREEDY + (Basic + Fix_A/Fix_B + LNS)* en boucle ===" << std::endl;
    std::cout << "=== Basic = solve_basic, Fix_A/B = lns_step_fix_A/B, LNS = lns_step_partial ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    auto [M, k, filename] = load_csv_matrix(csv_file, k_value);
    if (M.rows == 0) return result;
    
    auto total_start = std::chrono::high_resolution_clock::now();
    double timeout_ms = timeout_seconds * 1000.0;
    
    auto time_remaining_ms = [&]() -> double {
        if (timeout_ms <= 0) return std::numeric_limits<double>::max();
        auto now = std::chrono::high_resolution_clock::now();
        return timeout_ms - std::chrono::duration<double, std::milli>(now - total_start).count();
    };
    
    if (timeout_seconds > 0) {
        std::cout << "  Timeout: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "  Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    // Phase 1: GREEDY
    std::cout << "\n[Phase 1] Initialisation GREEDY..." << std::endl;
    BMFLocalSearch ls(k, M, 42);
    ls.initialize_greedy();
    ls.compute_all_counts();
    int init_errors = ls.count_errors();
    std::cout << "  [GREEDY] Erreurs: " << init_errors << std::endl;
    
    // Boucle (Basic + Fix_A/Fix_B + LNS)*
    int best_errors = init_errors;
    Matrix A_best = ls.A;
    Matrix B_best = ls.B;
    Matrix A_current = ls.A;
    Matrix B_current = ls.B;
    int total_loops = 0;
    int total_iterations = 0;
    std::mt19937 rng(42);
    int neighborhood_size = 10;
    int max_neighborhood = std::min(100, std::min(M.rows, M.cols));
    std::cout << "  Voisinage: " << neighborhood_size << " -> " << max_neighborhood 
              << " (par pas de 10, matrice " << M.rows << "x" << M.cols << ")" << std::endl;
    
    BMF bmf(k, M);
    
    while (best_errors > 0 && !is_interrupted() && time_remaining_ms() > 1000) {
        total_loops++;
        bool loop_improved = false;
        
        std::cout << "\n=== Boucle " << total_loops << " (erreurs: " << best_errors << ", voisinage=" << neighborhood_size << ") ===" << std::endl;
        
        // Phase A: Basic (solve_basic)
        std::cout << "[Basic] solve_basic jusqu'au minimum local..." << std::endl;
        BMFLocalSearch ls_basic(k, M, 42 + total_loops);
        ls_basic.A = A_current;
        ls_basic.B = B_current;
        ls_basic.compute_all_counts();
        
        LocalSearchResult basic_result = ls_basic.solve_basic(0, true);
        total_iterations += basic_result.iterations;
        
        if (basic_result.final_errors < best_errors) {
            best_errors = basic_result.final_errors;
            A_best = basic_result.A_solution;
            B_best = basic_result.B_solution;
            loop_improved = true;
            std::cout << "  Basic: " << best_errors << " erreurs" << std::endl;
        }
        
        A_current = basic_result.A_solution;
        B_current = basic_result.B_solution;
        
        if (best_errors == 0) break;
        if (is_interrupted() || time_remaining_ms() < 1000) break;
        
        // Phase B: Fix_A PROGRESSIF (colonne par colonne) - 70% top-k + 30% random
        std::vector<int> selected_cols = select_cols_topk_random(M, A_current, B_current, neighborhood_size, rng);
        std::cout << "[Fix_A] " << selected_cols.size() << " cols..." << std::endl;
        
        for (int col_j : selected_cols) {
            std::vector<int> single_col = {col_j};
            int new_cost = bmf.lns_step_fix_A(A_current, B_current, single_col);
            total_iterations++;
            if (new_cost >= 0 && new_cost < best_errors) {
                best_errors = new_cost;
                A_best = A_current;
                B_best = B_current;
                loop_improved = true;
                std::cout << "  \U0001F3AF Fix_A col " << col_j << ": " << best_errors << " erreurs" << std::endl;
                if (best_errors == 0) break;
            }
        }
        
        if (best_errors == 0) break;
        if (is_interrupted() || time_remaining_ms() < 1000) break;
        
        // Phase C: Fix_B PROGRESSIF (ligne par ligne) - 70% top-k + 30% random
        std::vector<int> selected_rows = select_rows_topk_random(M, A_current, B_current, neighborhood_size, rng);
        std::cout << "[Fix_B] " << selected_rows.size() << " rows..." << std::endl;
        
        for (int row_i : selected_rows) {
            std::vector<int> single_row = {row_i};
            int new_cost = bmf.lns_step_fix_B(A_current, B_current, single_row);
            total_iterations++;
            if (new_cost >= 0 && new_cost < best_errors) {
                best_errors = new_cost;
                A_best = A_current;
                B_best = B_current;
                loop_improved = true;
                std::cout << "  \U0001F3AF Fix_B row " << row_i << ": " << best_errors << " erreurs" << std::endl;
                if (best_errors == 0) break;
            }
        }
        
        if (best_errors == 0) break;
        if (is_interrupted() || time_remaining_ms() < 1000) break;
        
        // Phase D: LNS (lns_step_partial)
        std::cout << "[LNS] lns_step_partial (" << neighborhood_size << "x" << neighborhood_size << ")..." << std::endl;
        auto [rows_free, cols_free] = select_neighborhood_topk_random(M, A_current, B_current, neighborhood_size, rng);
        
        int new_cost = bmf.lns_step_partial(A_current, B_current, rows_free, cols_free);
        total_iterations++;
        
        if (new_cost >= 0 && new_cost < best_errors) {
            best_errors = new_cost;
            A_best = A_current;
            B_best = B_current;
            loop_improved = true;
            std::cout << "  LNS: " << best_errors << " erreurs" << std::endl;
        }
        
        A_current = A_best;
        B_current = B_best;
        
        // Détection de minimum local
        if (!loop_improved) {
            if (neighborhood_size >= max_neighborhood) {
                std::cout << "  Minimum local atteint à boucle " << total_loops 
                          << " avec " << best_errors << " erreurs (voisinage=" << max_neighborhood << " = max)" << std::endl;
                break;
            }
            int old_size = neighborhood_size;
            neighborhood_size = std::min(max_neighborhood, neighborhood_size + 10);
            std::cout << "  Pas d'amélioration, voisinage élargi: " 
                      << old_size << " -> " << neighborhood_size << std::endl;
        } else {
            neighborhood_size = 10;
        }
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(total_end - total_start).count();
    
    // Résultats
    std::cout << "\n=== RÉSULTATS OPTION 6 ===" << std::endl;
    std::cout << "Erreurs: " << init_errors << " -> " << best_errors << std::endl;
    std::cout << "Boucles: " << total_loops << std::endl;
    std::cout << "Itérations totales: " << total_iterations << std::endl;
    std::cout << "Temps total: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    
    result.k = k;
    result.init_errors = init_errors;
    result.final_errors = best_errors;
    result.iterations = total_iterations;
    result.total_time = total_time;
    result.success = true;
    
    int python_errors = verify_with_python(M, A_best, B_best, "OPT6_FULL_LOOP", filename);
    if (python_errors >= 0 && python_errors == best_errors) {
        std::cout << "Vérification OK" << std::endl;
    }
    
    return result;
}

/**
 * @brief Execute une option sur tous les fichiers CSV_FILES_WITH_K
 */
void run_batch_option(int option_num, double timeout_seconds) {
    std::cout << "\n=== MODE BATCH: Option " << option_num << " sur tous les fichiers ===" << std::endl;
    std::cout << "Timeout par fichier: " << timeout_seconds << " secondes" << std::endl;
    std::cout << "Fichier de résultats: " << RESULTS_OUTPUT_FILE << "\n" << std::endl;
    
    std::ofstream results_file(RESULTS_OUTPUT_FILE);
    if (!results_file.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir " << RESULTS_OUTPUT_FILE << std::endl;
        return;
    }
    
    results_file << "=== BATCH Option " << option_num << " ===\n\n";
    
    for (const auto& [filename, k] : CSV_FILES_WITH_K) {
        std::string full_path = CSV_BASE_DIR + filename;
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << ">>> " << filename << " (k=" << k << ")" << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
        CSVTestResult result;
        switch (option_num) {
            case 1: result = test_opt1_greedy_basic(full_path, k, timeout_seconds); break;
            case 2: result = test_opt2_greedy_fix_ab(full_path, k, timeout_seconds); break;
            case 3: result = test_opt3_greedy_lns_loop(full_path, k, timeout_seconds); break;
            case 4: result = test_opt4_greedy_basic_lns_loop(full_path, k, timeout_seconds); break;
            case 5: result = test_opt5_greedy_fix_lns_loop(full_path, k, timeout_seconds); break;
            case 6: result = test_opt6_greedy_full_loop(full_path, k, timeout_seconds); break;
            default: continue;
        }
        
        if (result.success) {
            results_file << filename << " k=" << k << " : "
                        << result.init_errors << " -> " << result.final_errors 
                        << " erreurs (" << std::fixed << std::setprecision(1) 
                        << result.total_time << " ms)\n";
            results_file.flush();
        }
        
        if (is_interrupted()) {
            std::cout << "\nBatch interrompu par l'utilisateur" << std::endl;
            break;
        }
    }
    
    results_file.close();
    std::cout << "\nRésultats sauvegardés dans: " << RESULTS_OUTPUT_FILE << std::endl;
}


// ==================== TEST LOOP: (LS + WLS + LNS-v3)* ====================

/**
 * @brief Test la factorisation BMF avec boucle externe (LS + WLS + LNS-v3)*
 * 
 * Exécute LS+WLS+LNS-v3 en boucle. Quand LNS-v3 stagne, on redémarre avec
 * la meilleure solution trouvée. Continue jusqu'à Ctrl+C, solution parfaite,
 * ou timeout global.
 * 
 * @param csv_file Chemin vers le fichier CSV
 * @param k_value Nombre de facteurs (0 = automatique)
 * @param timeout_seconds Timeout global en secondes (0 = pas de timeout, Ctrl+C pour arrêter)
 * @return CSVTestResult contenant les résultats du meilleur run
 */
CSVTestResult test_csv_loop(const std::string& csv_file, int k_value = 0, double timeout_seconds = 0) {
    CSVTestResult result;
    result.filename = std::filesystem::path(csv_file).filename().string();
    result.k = k_value;
    result.success = false;
    result.final_errors = INT_MAX;
    
    // Installer le gestionnaire de signal pour Ctrl+C
    std::signal(SIGINT, signal_handler);
    reset_interrupt();
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "=== BOUCLE (LS + WLS + LNS-v3)* - Redémarrage automatique ===" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::string filename = csv_file;
    
    // Si c'est juste un nom de fichier, ajouter le chemin data/
    if (filename.find('/') == std::string::npos) {
        filename = "../data/" + filename;
        if (filename.find(".csv") == std::string::npos) {
            filename += ".csv";
        }
    }
    
    std::cout << "Chargement de: " << filename << std::endl;
    
    Matrix M = CSVMatrixLoader::loadFromCSV(filename);
    DataGenerator::print_matrix_stats(M, "Matrice CSV");
    if (M.rows == 0 || M.cols == 0) {
        std::cerr << "Erreur: impossible de charger la matrice" << std::endl;
        return result;
    }
    
    int m = M.rows;
    int n = M.cols;
    
    // Déterminer k automatiquement si non spécifié
    int k = k_value;
    if (k <= 0) {
        k = std::max(5, std::min(50, static_cast<int>(std::log2(std::min(m, n)) * 2)));
        std::cout << "  k automatique: " << k << " (basé sur dimensions " << m << "x" << n << ")" << std::endl;
    } else {
        std::cout << "  k spécifié: " << k << std::endl;
    }
    
    if (timeout_seconds > 0) {
        std::cout << "  Timeout global: " << timeout_seconds << " secondes" << std::endl;
    } else {
        std::cout << "  Pas de timeout (Ctrl+C pour arrêter)" << std::endl;
    }
    
    // Timer global
    auto global_start = std::chrono::high_resolution_clock::now();
    double timeout_ms = timeout_seconds * 1000.0;
    
    auto time_elapsed_ms = [&]() -> double {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(now - global_start).count();
    };
    
    auto time_remaining_ms = [&]() -> double {
        if (timeout_ms <= 0) return std::numeric_limits<double>::max();
        return timeout_ms - time_elapsed_ms();
    };
    
    // Meilleure solution globale
    int best_errors = INT_MAX;
    Matrix A_best(0, 0), B_best(0, 0);
    int total_loops = 0;
    int total_lns_iterations = 0;
    double total_wls_time = 0;
    double total_lns_time = 0;
    
    // Solution courante (initialisée au premier passage)
    bool has_initial_solution = false;
    Matrix A_current(0, 0), B_current(0, 0);
    
    // ==================== BOUCLE PRINCIPALE ====================
    while (!is_interrupted() && time_remaining_ms() > 1000) {
        total_loops++;
        
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << ">>> BOUCLE " << total_loops << " (meilleur: " << best_errors 
                  << " erreurs, temps: " << std::fixed << std::setprecision(1) 
                  << time_elapsed_ms() / 1000.0 << "s)" << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
        // Varier le seed à chaque boucle
        //int seed = 42 + total_loops;
        int seed = 42;  // Multiplier pour plus de dispersion
        
        // === PHASE 1: GREEDY ===
        std::cout << "\n[Phase 1] Initialisation GREEDY (seed=" << seed << ")..." << std::endl;
        auto ls_start = std::chrono::high_resolution_clock::now();
        
        BMFLocalSearch ls(k, M, seed);
        
        if (has_initial_solution && best_errors < INT_MAX) {
            // Réutiliser la meilleure solution comme point de départ
            // Mais on perturbe légèrement pour explorer de nouvelles zones
            ls.A = A_best;
            ls.B = B_best;
            ls.compute_all_counts();
            
            // Perturbation légère
            /*std::uniform_int_distribution<> dist_m(0, m - 1);
            std::uniform_int_distribution<> dist_n(0, n - 1);
            std::uniform_int_distribution<> dist_k(0, k - 1);
            
            int perturb_count = std::max(1, best_errors / 10);  // ~10% des erreurs
            for (int p = 0; p < perturb_count; p++) {
                if (ls.rng() % 2 == 0) {
                    ls.flip_A(dist_m(ls.rng), dist_k(ls.rng));
                } else {
                    ls.flip_B(dist_k(ls.rng), dist_n(ls.rng));
                }
            }
            ls.compute_all_counts();
            std::cout << "  [RESTART] Réutilisation de la meilleure solution avec perturbation" << std::endl;*/
             std::cout << "  [RESTART] Réutilisation de la meilleure solution sans perturbation" << std::endl;
        } else {
            ls.initialize_greedy();
            ls.compute_all_counts();
        }
        
        int init_errors = ls.count_errors();
        std::cout << "  [GREEDY] Erreurs: " << init_errors << std::endl;
        
        // === PHASE 2: WLS ===
        int wls_errors = init_errors;
        if (init_errors > 0 && !is_interrupted() && time_remaining_ms() > 1000) {
            std::cout << "\n[Phase 2] WEIGHTED Local Search..." << std::endl;
            
            LocalSearchResult wls_result = ls.solve_weighted(
                100000,   // max_iterations
                0.2,        // penalty_increment
                30,         // max_stagnation
                true,       // verbose
                &g_interrupted
            );
            
            wls_errors = wls_result.final_errors;
            ls.A = wls_result.A_solution;
            ls.B = wls_result.B_solution;
            
            auto wls_end = std::chrono::high_resolution_clock::now();
            double wls_time = std::chrono::duration<double, std::milli>(wls_end - ls_start).count();
            total_wls_time += wls_time;
            
            std::cout << "  [WLS] " << init_errors << " -> " << wls_errors << " erreurs" << std::endl;
        }
        
        // Mettre à jour le meilleur si WLS a trouvé mieux
        if (wls_errors < best_errors) {
            best_errors = wls_errors;
            A_best = ls.A;
            B_best = ls.B;
            has_initial_solution = true;
            std::cout << "  NOUVEAU MEILLEUR après WLS: " << best_errors << " erreurs" << std::endl;
        }
        
        if (wls_errors == 0) {
            std::cout << "\nSolution parfaite trouvée par WLS!" << std::endl;
            break;
        }
        
        // === PHASE 3: LNS-v3 ===
        if (wls_errors > 0 && !is_interrupted() && time_remaining_ms() > 1000) {
            std::cout << "\n[Phase 3] LNS-v3 (Alternance A/B + Anti-stagnation)..." << std::endl;
            
            auto lns_start = std::chrono::high_resolution_clock::now();
            
            BMF bmf(k, M);
            bmf.set_initial_solution(ls.A, ls.B);
            
            // Paramètres LNS-v3 adaptés pour la boucle
            double neighborhood_ratio = 0.05;
            int max_iter_lns = 1000;
            int stagnation_threshold = 10;
            
            // Timeout par boucle LNS: temps FIXE par boucle (pas cumulatif!)
            // On alloue 60s max par phase LNS, ou le temps restant si moins
            double lns_timeout_per_loop = 60000.0;  // 60 secondes par boucle LNS
            double lns_timeout = std::min(lns_timeout_per_loop, time_remaining_ms() - 1000.0);
            
            auto lns_result = bmf.solve_lns_v3(
                max_iter_lns,
                neighborhood_ratio,
                seed,
                stagnation_threshold,
                true,
                lns_timeout,
                &g_interrupted
            );
            
            auto lns_end = std::chrono::high_resolution_clock::now();
            double lns_time = std::chrono::duration<double, std::milli>(lns_end - lns_start).count();
            total_lns_time += lns_time;
            total_lns_iterations += lns_result.iterations;
            
            int lns_errors = lns_result.final_errors;
            
            std::cout << "  [LNS-v3] " << wls_errors << " -> " << lns_errors 
                      << " erreurs (" << lns_result.iterations << " iters, "
                      << lns_result.stop_reason << ")" << std::endl;
            
            // Mettre à jour le meilleur si LNS a trouvé mieux
            if (lns_errors < best_errors) {
                best_errors = lns_errors;
                A_best = lns_result.A_solution;
                B_best = lns_result.B_solution;
                has_initial_solution = true;
                std::cout << "  NOUVEAU MEILLEUR après LNS-v3: " << best_errors << " erreurs" << std::endl;
            }
            
            if (lns_errors == 0) {
                std::cout << "\nSolution parfaite trouvée par LNS-v3!" << std::endl;
                break;
            }
        }
        
        // Afficher le résumé de la boucle
        std::cout << "\n--- Fin boucle " << total_loops << ": meilleur = " << best_errors 
                  << " erreurs ---" << std::endl;
    }
    
    // ==================== RÉSUMÉ FINAL ====================
    auto global_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(global_end - global_start).count();
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== RÉSUMÉ BOUCLE (LS + WLS + LNS-v3)* ===" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "Fichier: " << CSVMatrixLoader::getBaseName(filename) << std::endl;
    std::cout << "Dimensions: " << m << "x" << n << ", k=" << k << std::endl;
    std::cout << "Nombre de boucles: " << total_loops << std::endl;
    std::cout << "Itérations LNS-v3 totales: " << total_lns_iterations << std::endl;
    std::cout << "Temps WLS total: " << std::fixed << std::setprecision(1) << total_wls_time << " ms" << std::endl;
    std::cout << "Temps LNS-v3 total: " << std::fixed << std::setprecision(1) << total_lns_time << " ms" << std::endl;
    std::cout << "Temps total: " << std::fixed << std::setprecision(1) << total_time << " ms" << std::endl;
    std::cout << "\nMEILLEUR RÉSULTAT: " << best_errors << " erreurs" << std::endl;
    
    // Raison d'arrêt
    std::string stop_reason = (best_errors == 0) ? "zero_errors" :
                              is_interrupted() ? "interrupted" :
                              (time_remaining_ms() <= 0) ? "timeout" : "unknown";
    std::string stop_label = (stop_reason == "zero_errors") ? "[OK]" :
                             (stop_reason == "interrupted") ? "[INTERROMPU]" :
                             (stop_reason == "timeout") ? "[TIMEOUT]" : "[INCONNU]";
    std::cout << "Arret: " << stop_label << " " << stop_reason << std::endl;
    
    // Remplir le résultat
    result.k = k;
    result.init_errors = 0;  // N/A pour loop
    result.after_ls_errors = 0;  // N/A pour loop
    result.final_errors = best_errors;
    result.iterations = total_lns_iterations;
    result.ls_time = total_wls_time;
    result.lns_time = total_lns_time;
    result.total_time = total_time;
    result.success = (best_errors == 0);
    
    // ==================== VÉRIFICATION PYTHON ====================
    if (best_errors < INT_MAX) {
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "=== VÉRIFICATION PYTHON ===" << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
        int python_errors = verify_with_python(M, A_best, B_best, "LOOP", filename);
        
        if (python_errors >= 0) {
            if (python_errors == best_errors) {
                std::cout << "Vérification OK: C++ (" << best_errors << ") == Python (" << python_errors << ")" << std::endl;
            } else {
                std::cout << "Différence: C++ (" << best_errors << ") != Python (" << python_errors << ")" << std::endl;
            }
        }
    }
    
    return result;
}


// ==================== MAIN ====================

void print_usage() {
    std::cout << "Usage: ./rbac_main [option]" << std::endl;
    std::cout << "\n=== OPTIONS PRINCIPALES ===" << std::endl;
    std::cout << "  --csv <file> <k> [timeout]       Test GREEDY + WLS + LNS-v3 sur un fichier CSV" << std::endl;
    std::cout << "  --csv [timeout]                  Test batch sur tous les fichiers prédéfinis" << std::endl;
    std::cout << "  --alt <file> <k> [timeout]       Test avec alternance (WLS <-> LNS-v3)*" << std::endl;
    std::cout << "  --multi <file> <k> <N> [timeout] Multi-start: N exécutions avec seeds différents" << std::endl;
    std::cout << "  --loop <file> <k> [timeout]      Boucle (LS+WLS+LNS-v3)* avec redémarrage auto" << std::endl;
    std::cout << "\n=== OPTIONS DE TEST MODULAIRES ===" << std::endl;
    std::cout << "  --opt1 <file> <k> [timeout]      GREEDY + solve_basic (jusqu'au min local)" << std::endl;
    std::cout << "  --opt2 <file> <k> [timeout]      GREEDY + Fix_A/Fix_B (lns_step_fix_A + lns_step_fix_B)" << std::endl;
    std::cout << "  --opt3 <file> <k> [timeout]      GREEDY + (LNS)* en boucle (lns_step_partial)" << std::endl;
    std::cout << "  --opt4 <file> <k> [timeout]      GREEDY + (solve_basic + lns_step_partial)* en boucle" << std::endl;
    std::cout << "  --opt5 <file> <k> [timeout]      GREEDY + (Fix_A/Fix_B + lns_step_partial)* en boucle" << std::endl;
    std::cout << "  --opt6 <file> <k> [timeout]      GREEDY + (solve_basic + Fix_A/Fix_B + lns_step_partial)*" << std::endl;
    std::cout << "\n=== AUTRES ===" << std::endl;
    std::cout << "  --lswls <file> <k>               Test GREEDY + WLS seulement (sans LNS-v3)" << std::endl;
    std::cout << "  --lslns <file> <k> [timeout]     Test GREEDY + LNS-v3 directement (sans WLS)" << std::endl;
    std::cout << "  --help                           Afficher cette aide" << std::endl;
    std::cout << "\n=== NOTES ===" << std::endl;
    std::cout << "  - Si <file> n'est pas spécifié, traite tous les fichiers de CSV_FILES_WITH_K" << std::endl;
    std::cout << "  - timeout = 0 signifie pas de limite (Ctrl+C pour arrêter)" << std::endl;
    std::cout << "  - Les résultats batch sont sauvegardés dans: " << RESULTS_OUTPUT_FILE << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "================================================" << std::endl;
    std::cout << "  BMF Solver (Boolean Matrix Factorization)" << std::endl;
    std::cout << "  Algorithmes: GREEDY + WLS + LNS-v3" << std::endl;
    std::cout << "================================================" << std::endl;
    
    if (argc > 1) {
        std::string arg = argv[1];
        
        if (arg == "--csv") {
            if (argc >= 3) {
                std::string arg2 = argv[2];
                // Vérifier si arg2 est un nombre (timeout) ou un fichier
                bool is_number = !arg2.empty() && std::all_of(arg2.begin(), arg2.end(), [](char c) {
                    return std::isdigit(c) || c == '.';
                });
                
                if (is_number) {
                    // Mode batch: --csv [timeout]
                    double timeout = std::stod(arg2);
                    std::cout << "\n=== MODE BATCH: Exécution sur les fichiers prédéfinis ===" << std::endl;
                    std::cout << "Timeout par fichier: " << timeout << " secondes" << std::endl;
                    std::cout << "Fichier de résultats: " << RESULTS_OUTPUT_FILE << "\n" << std::endl;
                    
                    std::ofstream results_file(RESULTS_OUTPUT_FILE);
                    if (!results_file.is_open()) {
                        std::cerr << "Erreur: Impossible d'ouvrir " << RESULTS_OUTPUT_FILE << std::endl;
                        return 1;
                    }
                    
                    for (const auto& [filename, k] : CSV_FILES_WITH_K) {
                        std::string full_path = CSV_BASE_DIR + filename;
                        std::cout << "\n" << std::string(70, '=') << std::endl;
                        std::cout << ">>> " << filename << " (k=" << k << ")" << std::endl;
                        std::cout << std::string(70, '=') << std::endl;
                        
                        CSVTestResult result = test_csv_partial(full_path, k, timeout);
                        
                        if (result.success) {
                            results_file << result.to_string() << "\n";
                            results_file.flush();
                        }
                        
                        if (is_interrupted()) {
                            std::cout << "\nBatch interrompu par l'utilisateur" << std::endl;
                            break;
                        }
                    }
                    
                    results_file.close();
                    std::cout << "\n Résultats sauvegardés dans: " << RESULTS_OUTPUT_FILE << std::endl;
                } else {
                    // Mode fichier unique: --csv <fichier> [k] [timeout]
                    std::string csv_file = arg2;
                    int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                    double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                    test_csv_partial(csv_file, k, timeout);
                }
            } else {
                // Mode batch sans timeout: --csv (seul)
                std::cout << "\n=== MODE BATCH: Exécution sur les fichiers prédéfinis ===" << std::endl;
                std::cout << "Sans timeout (Ctrl+C pour arrêter)" << std::endl;
                std::cout << "Fichier de résultats: " << RESULTS_OUTPUT_FILE << "\n" << std::endl;
                
                std::ofstream results_file(RESULTS_OUTPUT_FILE);
                if (!results_file.is_open()) {
                    std::cerr << "Erreur: Impossible d'ouvrir " << RESULTS_OUTPUT_FILE << std::endl;
                    return 1;
                }
                
                for (const auto& [filename, k] : CSV_FILES_WITH_K) {
                    std::string full_path = CSV_BASE_DIR + filename;
                    std::cout << "\n" << std::string(70, '=') << std::endl;
                    std::cout << ">>> " << filename << " (k=" << k << ")" << std::endl;
                    std::cout << std::string(70, '=') << std::endl;
                    
                    CSVTestResult result = test_csv_partial(full_path, k, 0);
                    
                    if (result.success) {
                        results_file << result.to_string() << "\n";
                        results_file.flush();
                    }
                    
                    if (is_interrupted()) {
                        std::cout << "\nBatch interrompu par l'utilisateur" << std::endl;
                        break;
                    }
                }
                
                results_file.close();
                std::cout << "\nRésultats sauvegardés dans: " << RESULTS_OUTPUT_FILE << std::endl;
            }
        } else if (arg == "--alt") {
            // Mode alternance WLS <-> LNS-v3
            if (argc >= 3) {
                std::string csv_file = argv[2];
                int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                test_csv_alternating(csv_file, k, timeout);
            } else {
                std::cout << "Erreur: --alt nécessite au moins un fichier CSV" << std::endl;
                std::cout << "Usage: --alt <file> <k> [timeout]" << std::endl;
            }
        } else if (arg == "--multi") {
            // Mode multi-start: N exécutions avec seeds différents
            if (argc >= 5) {
                std::string csv_file = argv[2];
                int k = std::stoi(argv[3]);
                int num_starts = std::stoi(argv[4]);
                double timeout = (argc >= 6) ? std::stod(argv[5]) : 0;
                test_csv_multistart(csv_file, k, num_starts, timeout);
            } else {
                std::cout << "Erreur: --multi nécessite fichier, k et N" << std::endl;
                std::cout << "Usage: --multi <file> <k> <N> [timeout]" << std::endl;
                std::cout << "  N = nombre d'exécutions avec seeds différents" << std::endl;
            }
        } else if (arg == "--lswls") {
            // Mode LS + WLS seulement (sans LNS-v3)
            if (argc >= 3) {
                std::string csv_file = argv[2];
                int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                test_csv_lswls(csv_file, k);
            } else {
                std::cout << "Erreur: --lswls nécessite au moins un fichier CSV" << std::endl;
                std::cout << "Usage: --lswls <file> <k>" << std::endl;
            }
        } else if (arg == "--lslns") {
            // Mode LS + LNS-v3 directement (sans WLS)
            if (argc >= 3) {
                std::string csv_file = argv[2];
                int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                test_csv_lslns(csv_file, k, timeout);
            } else {
                std::cout << "Erreur: --lslns nécessite au moins un fichier CSV" << std::endl;
                std::cout << "Usage: --lslns <file> <k> [timeout]" << std::endl;
            }
        } else if (arg == "--loop") {
            // Mode boucle (LS + WLS + LNS-v3)* avec redémarrage automatique
            if (argc >= 3) {
                std::string csv_file = argv[2];
                int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                test_csv_loop(csv_file, k, timeout);
            } else {
                std::cout << "Erreur: --loop nécessite au moins un fichier CSV" << std::endl;
                std::cout << "Usage: --loop <file> <k> [timeout]" << std::endl;
            }
        } else if (arg == "--opt1") {
            // Option 1: GREEDY + Local Search Basic
            if (argc >= 3) {
                std::string arg2 = argv[2];
                bool is_number = !arg2.empty() && std::all_of(arg2.begin(), arg2.end(), [](char c) {
                    return std::isdigit(c) || c == '.';
                });
                if (is_number) {
                    // Mode batch
                    run_batch_option(1, std::stod(arg2));
                } else {
                    std::string csv_file = arg2;
                    int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                    double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                    test_opt1_greedy_basic(csv_file, k, timeout);
                }
            } else {
                run_batch_option(1, 0);  // Batch sans timeout
            }
        } else if (arg == "--opt2") {
            // Option 2: GREEDY + WLS
            if (argc >= 3) {
                std::string arg2 = argv[2];
                bool is_number = !arg2.empty() && std::all_of(arg2.begin(), arg2.end(), [](char c) {
                    return std::isdigit(c) || c == '.';
                });
                if (is_number) {
                    run_batch_option(2, std::stod(arg2));
                } else {
                    std::string csv_file = arg2;
                    int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                    double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                    test_opt2_greedy_fix_ab(csv_file, k, timeout);
                }
            } else {
                run_batch_option(2, 0);
            }
        } else if (arg == "--opt3") {
            // Option 3: GREEDY + (LNS)* en boucle
            if (argc >= 3) {
                std::string arg2 = argv[2];
                bool is_number = !arg2.empty() && std::all_of(arg2.begin(), arg2.end(), [](char c) {
                    return std::isdigit(c) || c == '.';
                });
                if (is_number) {
                    run_batch_option(3, std::stod(arg2));
                } else {
                    std::string csv_file = arg2;
                    int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                    double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                    test_opt3_greedy_lns_loop(csv_file, k, timeout);
                }
            } else {
                run_batch_option(3, 0);
            }
        } else if (arg == "--opt4") {
            // Option 4: GREEDY + (Basic + LNS)* en boucle
            if (argc >= 3) {
                std::string arg2 = argv[2];
                bool is_number = !arg2.empty() && std::all_of(arg2.begin(), arg2.end(), [](char c) {
                    return std::isdigit(c) || c == '.';
                });
                if (is_number) {
                    run_batch_option(4, std::stod(arg2));
                } else {
                    std::string csv_file = arg2;
                    int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                    double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                    test_opt4_greedy_basic_lns_loop(csv_file, k, timeout);
                }
            } else {
                run_batch_option(4, 0);
            }
        } else if (arg == "--opt5") {
            // Option 5: GREEDY + (WLS + LNS)* en boucle
            if (argc >= 3) {
                std::string arg2 = argv[2];
                bool is_number = !arg2.empty() && std::all_of(arg2.begin(), arg2.end(), [](char c) {
                    return std::isdigit(c) || c == '.';
                });
                if (is_number) {
                    run_batch_option(5, std::stod(arg2));
                } else {
                    std::string csv_file = arg2;
                    int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                    double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                    test_opt5_greedy_fix_lns_loop(csv_file, k, timeout);
                }
            } else {
                run_batch_option(5, 0);
            }
        } else if (arg == "--opt6") {
            // Option 6: GREEDY + (Basic + WLS + LNS)* en boucle
            if (argc >= 3) {
                std::string arg2 = argv[2];
                bool is_number = !arg2.empty() && std::all_of(arg2.begin(), arg2.end(), [](char c) {
                    return std::isdigit(c) || c == '.';
                });
                if (is_number) {
                    run_batch_option(6, std::stod(arg2));
                } else {
                    std::string csv_file = arg2;
                    int k = (argc >= 4) ? std::stoi(argv[3]) : 0;
                    double timeout = (argc >= 5) ? std::stod(argv[4]) : 0;
                    test_opt6_greedy_full_loop(csv_file, k, timeout);
                }
            } else {
                run_batch_option(6, 0);
            }
        } else if (arg == "--help") {
            print_usage();
        } else {
            std::cout << "Option inconnue: " << arg << std::endl;
            print_usage();
        }
    } else {
        print_usage();
    }
    
    return 0;
}
