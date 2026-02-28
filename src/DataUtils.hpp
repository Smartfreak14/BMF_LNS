#pragma once
#include "Matrix.hpp"
#include <random>
#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <map>
#include <vector>
#include <optional>
#include <regex>

// ==================== JSON SIMPLE PARSER/WRITER ====================
// Note: Pour un projet plus robuste, utiliser nlohmann/json

namespace json {

/**
 * @brief Écrit une matrice au format JSON
 */
inline std::string matrix_to_json(const Matrix& M) {
    std::ostringstream oss;
    oss << "[";
    for (int i = 0; i < M.rows; i++) {
        oss << "[";
        for (int j = 0; j < M.cols; j++) {
            oss << M(i, j);
            if (j < M.cols - 1) oss << ",";
        }
        oss << "]";
        if (i < M.rows - 1) oss << ",";
    }
    oss << "]";
    return oss.str();
}

/**
 * @brief Parse un tableau JSON simple [[1,0],[0,1]] en matrice
 */
inline Matrix json_to_matrix(const std::string& json_str) {
    std::vector<std::vector<int>> data;
    std::vector<int> current_row;
    std::string num_buffer;
    int depth = 0;
    
    for (char c : json_str) {
        if (c == '[') {
            depth++;
            if (depth == 2) {
                current_row.clear();
            }
        } else if (c == ']') {
            if (!num_buffer.empty()) {
                current_row.push_back(std::stoi(num_buffer));
                num_buffer.clear();
            }
            if (depth == 2 && !current_row.empty()) {
                data.push_back(current_row);
            }
            depth--;
        } else if (c == ',') {
            if (!num_buffer.empty()) {
                current_row.push_back(std::stoi(num_buffer));
                num_buffer.clear();
            }
        } else if (c >= '0' && c <= '9') {
            num_buffer += c;
        }
    }
    
    if (data.empty()) return Matrix(0, 0);
    return Matrix(data);
}

} // namespace json

/**
 * @class DataGenerator
 * @brief Génère des configurations RBAC aléatoires pour les tests
 */
class DataGenerator {
public:
    std::mt19937 rng;
    
    DataGenerator(unsigned seed = std::random_device{}()) : rng(seed) {}
    
    /**
     * @brief Génère une matrice aléatoire avec une densité donnée
     */
    Matrix random_matrix(int rows, int cols, double density = 0.5) {
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
    
    /**
     * @brief Génère une configuration RBAC (A, B)
     * @param m Nombre d'utilisateurs
     * @param n_objects Nombre d'objets
     * @param k Nombre de rôles
     * @param density Densité des 1 dans les matrices
     */
    std::pair<Matrix, Matrix> generate_rbac_config(int m, int n_objects, int k, double density = 0.3) {
        Matrix A = random_matrix(m, k, density);
        Matrix B = random_matrix(k, n_objects * 2, density);  // 2 colonnes par objet (read, write)
        return {A, B};
    }
    
    /**
     * @brief Génère des permissions aléatoires pour un nouveau sujet
     */
    std::vector<int> random_subject_permissions(int n_objects, double density = 0.3) {
        std::vector<int> perms(n_objects * 2, 0);
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        for (int i = 0; i < n_objects * 2; i++) {
            if (dis(rng) < density) {
                perms[i] = 1;
            }
        }
        return perms;
    }
    
    /**
     * @brief Génère des permissions aléatoires pour un nouvel objet
     */
    std::vector<std::vector<int>> random_object_permissions(int m, double density = 0.3) {
        std::vector<std::vector<int>> perms(m, {0, 0});
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        for (int i = 0; i < m; i++) {
            if (dis(rng) < density) perms[i][0] = 1;  // read
            if (dis(rng) < density) perms[i][1] = 1;  // write
        }
        return perms;
    }
};

/**
 * @class ConfigurationManager
 * @brief Sauvegarde et charge des configurations RBAC au format CSV/JSON
 */
class ConfigurationManager {
public:
    
    /**
     * @brief Sauvegarde une configuration RBAC en CSV
     */
    static bool save_csv(const Matrix& A, const Matrix& B, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir " << filename << std::endl;
            return false;
        }
        
        // Header
        file << "type,row,";
        for (int j = 0; j < std::max(A.cols, B.cols); j++) {
            file << "col" << j;
            if (j < std::max(A.cols, B.cols) - 1) file << ",";
        }
        file << "\n";
        
        // Matrix A
        for (int i = 0; i < A.rows; i++) {
            file << "A," << i << ",";
            for (int j = 0; j < A.cols; j++) {
                file << A(i, j);
                if (j < A.cols - 1) file << ",";
            }
            file << "\n";
        }
        
        // Matrix B
        for (int i = 0; i < B.rows; i++) {
            file << "B," << i << ",";
            for (int j = 0; j < B.cols; j++) {
                file << B(i, j);
                if (j < B.cols - 1) file << ",";
            }
            file << "\n";
        }
        
        file.close();
        std::cout << "Configuration sauvegardée dans " << filename << std::endl;
        return true;
    }
    
    /**
     * @brief Charge une configuration RBAC depuis un CSV
     */
    static std::pair<Matrix, Matrix> load_csv(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir " << filename << std::endl;
            return {Matrix(0, 0), Matrix(0, 0)};
        }
        
        std::vector<std::vector<int>> A_data, B_data;
        std::string line;
        
        // Skip header
        std::getline(file, line);
        
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string type, row_str, val;
            
            std::getline(ss, type, ',');
            std::getline(ss, row_str, ',');
            
            std::vector<int> row;
            while (std::getline(ss, val, ',')) {
                if (!val.empty()) {
                    row.push_back(std::stoi(val));
                }
            }
            
            if (type == "A") {
                A_data.push_back(row);
            } else if (type == "B") {
                B_data.push_back(row);
            }
        }
        
        file.close();
        
        return {Matrix(A_data), Matrix(B_data)};
    }
    
    /**
     * @brief Génère un nom de fichier avec timestamp
     */
    static std::string generate_filename(int m, int n_obj, int k, const std::string& prefix = "rbac_config") {
        auto now = std::time(nullptr);
        auto tm = *std::localtime(&now);
        
        std::ostringstream oss;
        oss << prefix << "_m" << m << "_n" << n_obj << "_k" << k << "_";
        oss << std::put_time(&tm, "%Y%m%d_%H%M%S") << ".csv";
        
        return oss.str();
    }
};

/**
 * @class ResultsRecorder
 * @brief Enregistre les résultats d'expériences
 */
class ResultsRecorder {
private:
    std::vector<std::vector<std::string>> records;
    std::vector<std::string> headers;
    
public:
    ResultsRecorder(const std::vector<std::string>& cols) : headers(cols) {}
    
    void add_record(const std::vector<std::string>& record) {
        records.push_back(record);
    }
    
    void save(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir " << filename << std::endl;
            return;
        }
        
        // Header
        for (size_t i = 0; i < headers.size(); i++) {
            file << headers[i];
            if (i < headers.size() - 1) file << ",";
        }
        file << "\n";
        
        // Data
        for (const auto& record : records) {
            for (size_t i = 0; i < record.size(); i++) {
                file << record[i];
                if (i < record.size() - 1) file << ",";
            }
            file << "\n";
        }
        
        file.close();
        std::cout << "Résultats sauvegardés dans " << filename << std::endl;
    }
    
    void print() {
        // Header
        for (const auto& h : headers) {
            std::cout << std::setw(15) << h;
        }
        std::cout << "\n" << std::string(headers.size() * 15, '-') << "\n";
        
        // Data
        for (const auto& record : records) {
            for (const auto& val : record) {
                std::cout << std::setw(15) << val;
            }
            std::cout << "\n";
        }
    }
};

// ==================== CONFIGURATION JSON MANAGER ====================

/**
 * @struct ExperimentParams
 * @brief Paramètres d'une expérience
 */
struct ExperimentParams {
    std::vector<int> subject_permissions;
    std::vector<std::vector<int>> object_permissions;
    int subject_id;
    int object_id;
    std::string operation;  // "read" ou "write"
};

/**
 * @struct RBACConfiguration
 * @brief Une configuration RBAC complète
 */
struct RBACConfiguration {
    int m, n_objects, k;
    std::string timestamp;
    Matrix A, B, M;
    std::map<int, ExperimentParams> experiments;
    
    RBACConfiguration() : m(0), n_objects(0), k(0), A(0,0), B(0,0), M(0,0) {}
    RBACConfiguration(int m_, int n_, int k_, const Matrix& A_, const Matrix& B_)
        : m(m_), n_objects(n_), k(k_), A(A_), B(B_), M(A_.multiply(B_)) {}
};

/**
 * @class JSONConfigurationManager
 * @brief Gère les configurations RBAC au format JSON (comme le code Python)
 */
class JSONConfigurationManager {
private:
    std::map<std::string, RBACConfiguration> configurations;
    std::string config_file;
    
public:
    JSONConfigurationManager(const std::string& filename = "rbac_all_configurations.json")
        : config_file(filename) {}
    
    /**
     * @brief Génère une clé de configuration
     */
    static std::string config_key(int m, int n_objects, int k) {
        return "m" + std::to_string(m) + "_n" + std::to_string(n_objects) + "_k" + std::to_string(k);
    }
    
    /**
     * @brief Charge les configurations depuis le fichier JSON
     */
    bool load() {
        std::ifstream file(config_file);
        if (!file.is_open()) {
            std::cout << "Fichier " << config_file << " non trouvé, sera créé." << std::endl;
            return false;
        }
        
        // Lire le contenu complet
        std::stringstream buffer;
        buffer << file.rdbuf();
        std::string content = buffer.str();
        file.close();
        
        // Parser le JSON (simplifiée - pour un vrai projet, utiliser nlohmann/json)
        // Cette implémentation simple cherche les patterns "configurations" et parse les matrices
        
        size_t configs_pos = content.find("\"configurations\"");
        if (configs_pos == std::string::npos) {
            std::cerr << "Format JSON invalide" << std::endl;
            return false;
        }
        
        // Parser chaque configuration m{m}_n{n}_k{k}
        std::regex config_regex("\"(m\\d+_n\\d+_k\\d+)\"\\s*:\\s*\\{");
        std::smatch match;
        std::string search_content = content.substr(configs_pos);
        
        while (std::regex_search(search_content, match, config_regex)) {
            std::string key = match[1].str();
            
            // Extraire m, n, k du nom
            int m, n, k_val;
            if (sscanf(key.c_str(), "m%d_n%d_k%d", &m, &n, &k_val) == 3) {
                // Trouver les matrices A et B
                size_t a_pos = search_content.find("\"A\":", match.position());
                size_t b_pos = search_content.find("\"B\":", match.position());
                
                if (a_pos != std::string::npos && b_pos != std::string::npos) {
                    // Extraire la matrice A
                    size_t a_start = search_content.find('[', a_pos);
                    int depth = 0;
                    size_t a_end = a_start;
                    for (size_t i = a_start; i < search_content.size(); i++) {
                        if (search_content[i] == '[') depth++;
                        else if (search_content[i] == ']') {
                            depth--;
                            if (depth == 0) { a_end = i + 1; break; }
                        }
                    }
                    std::string a_json = search_content.substr(a_start, a_end - a_start);
                    
                    // Extraire la matrice B
                    size_t b_start = search_content.find('[', b_pos);
                    depth = 0;
                    size_t b_end = b_start;
                    for (size_t i = b_start; i < search_content.size(); i++) {
                        if (search_content[i] == '[') depth++;
                        else if (search_content[i] == ']') {
                            depth--;
                            if (depth == 0) { b_end = i + 1; break; }
                        }
                    }
                    std::string b_json = search_content.substr(b_start, b_end - b_start);
                    
                    Matrix A = json::json_to_matrix(a_json);
                    Matrix B = json::json_to_matrix(b_json);
                    
                    if (A.rows > 0 && B.rows > 0) {
                        configurations[key] = RBACConfiguration(m, n, k_val, A, B);
                    }
                }
            }
            
            search_content = match.suffix().str();
        }
        
        std::cout << "" << configurations.size() << " configurations chargées depuis " << config_file << std::endl;
        return true;
    }
    
    /**
     * @brief Sauvegarde toutes les configurations dans le fichier JSON
     */
    bool save() {
        std::ofstream file(config_file);
        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir " << config_file << std::endl;
            return false;
        }
        
        auto now = std::time(nullptr);
        auto tm = *std::localtime(&now);
        
        file << "{\n";
        file << "  \"metadata\": {\n";
        file << "    \"created\": \"" << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S") << "\",\n";
        file << "    \"total_configs\": " << configurations.size() << "\n";
        file << "  },\n";
        file << "  \"configurations\": {\n";
        
        size_t count = 0;
        for (const auto& [key, config] : configurations) {
            file << "    \"" << key << "\": {\n";
            file << "      \"metadata\": {\n";
            file << "        \"m\": " << config.m << ",\n";
            file << "        \"n_objects\": " << config.n_objects << ",\n";
            file << "        \"k\": " << config.k << ",\n";
            file << "        \"satisfies_constraint\": true\n";
            file << "      },\n";
            file << "      \"matrices\": {\n";
            file << "        \"A\": " << json::matrix_to_json(config.A) << ",\n";
            file << "        \"B\": " << json::matrix_to_json(config.B) << ",\n";
            file << "        \"M\": " << json::matrix_to_json(config.M) << "\n";
            file << "      },\n";
            file << "      \"experiments\": {\n";
            
            size_t exp_count = 0;
            for (const auto& [exp_id, exp_params] : config.experiments) {
                file << "        \"exp_" << exp_id << "\": {\n";
                file << "          \"subject_permissions\": [";
                for (size_t i = 0; i < exp_params.subject_permissions.size(); i++) {
                    file << exp_params.subject_permissions[i];
                    if (i < exp_params.subject_permissions.size() - 1) file << ",";
                }
                file << "],\n";
                file << "          \"object_permissions\": [";
                for (size_t i = 0; i < exp_params.object_permissions.size(); i++) {
                    file << "[" << exp_params.object_permissions[i][0] << "," 
                         << exp_params.object_permissions[i][1] << "]";
                    if (i < exp_params.object_permissions.size() - 1) file << ",";
                }
                file << "],\n";
                file << "          \"subject_id\": " << exp_params.subject_id << ",\n";
                file << "          \"object_id\": " << exp_params.object_id << ",\n";
                file << "          \"operation\": \"" << exp_params.operation << "\"\n";
                file << "        }";
                if (++exp_count < config.experiments.size()) file << ",";
                file << "\n";
            }
            
            file << "      }\n";
            file << "    }";
            if (++count < configurations.size()) file << ",";
            file << "\n";
        }
        
        file << "  }\n";
        file << "}\n";
        
        file.close();
        std::cout << "" << configurations.size() << " configurations sauvegardées dans " << config_file << std::endl;
        return true;
    }
    
    /**
     * @brief Ajoute une configuration
     */
    void add_configuration(int m, int n_objects, int k, const Matrix& A, const Matrix& B) {
        std::string key = config_key(m, n_objects, k);
        configurations[key] = RBACConfiguration(m, n_objects, k, A, B);
    }
    
    /**
     * @brief Récupère une configuration
     */
    std::optional<RBACConfiguration> get_configuration(int m, int n_objects, int k) const {
        std::string key = config_key(m, n_objects, k);
        auto it = configurations.find(key);
        if (it != configurations.end()) {
            return it->second;
        }
        return std::nullopt;
    }
    
    /**
     * @brief Ajoute des paramètres d'expérience à une configuration
     */
    void add_experiment(int m, int n_objects, int k, int exp_id, const ExperimentParams& params) {
        std::string key = config_key(m, n_objects, k);
        auto it = configurations.find(key);
        if (it != configurations.end()) {
            it->second.experiments[exp_id] = params;
        }
    }
    
    /**
     * @brief Génère les paramètres d'expérience aléatoires
     */
    ExperimentParams generate_experiment_params(int m, int n_objects, unsigned seed) {
        std::mt19937 rng(seed);
        std::uniform_int_distribution<> perm_dist(0, 1);
        std::uniform_int_distribution<> subj_dist(0, m - 1);
        std::uniform_int_distribution<> obj_dist(0, n_objects - 1);
        std::uniform_int_distribution<> op_dist(0, 1);
        
        ExperimentParams params;
        
        // Subject permissions (2 * n_objects)
        params.subject_permissions.resize(n_objects * 2);
        for (int i = 0; i < n_objects * 2; i++) {
            params.subject_permissions[i] = perm_dist(rng);
        }
        
        // Object permissions (m utilisateurs)
        params.object_permissions.resize(m);
        for (int i = 0; i < m; i++) {
            params.object_permissions[i] = {perm_dist(rng), perm_dist(rng)};
        }
        
        params.subject_id = subj_dist(rng);
        params.object_id = obj_dist(rng);
        params.operation = (op_dist(rng) == 0) ? "read" : "write";
        
        return params;
    }
    
    /**
     * @brief Vérifie si une configuration existe
     */
    bool has_configuration(int m, int n_objects, int k) const {
        return configurations.find(config_key(m, n_objects, k)) != configurations.end();
    }
    
    /**
     * @brief Retourne le nombre de configurations
     */
    size_t size() const { return configurations.size(); }
    
    /**
     * @brief Liste les configurations
     */
    void list() const {
        std::cout << "Configurations disponibles (" << configurations.size() << "):\n";
        for (const auto& [key, config] : configurations) {
            std::cout << "   - " << key << " (m=" << config.m 
                      << ", n=" << config.n_objects 
                      << ", k=" << config.k << ")\n";
        }
    }
};

// ==================== RESULTS CSV MANAGER ====================

/**
 * @class ExperimentResultsRecorder
 * @brief Enregistre les résultats d'expériences au format CSV (comme le code Python)
 */
class ExperimentResultsRecorder {
public:
    struct ExperimentResult {
        int m, n_objects, k, exp_id;
        // Add-Subject
        double user_maxsat_time;
        int user_modif, user_a, user_b, user_success;
        // Add-Object
        double resource_maxsat_time;
        int resource_modif, resource_a, resource_b, resource_success;
        // Modify-Rule
        double modify_maxsat_time;
        int modify_modif, modify_a, modify_b, modify_success;
    };
    
private:
    std::vector<ExperimentResult> results;
    std::string output_file;
    
public:
    ExperimentResultsRecorder(const std::string& filename = "rbac_experiments_results.csv")
        : output_file(filename) {}
    
    void add_result(const ExperimentResult& result) {
        results.push_back(result);
    }
    
    bool save() {
        std::ofstream file(output_file);
        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir " << output_file << std::endl;
            return false;
        }
        
        // Header
        file << "m,n_objects,k,exp_id,"
             << "user_maxsat_time,user_modif,user_a,user_b,user_success,"
             << "resource_maxsat_time,resource_modif,resource_a,resource_b,resource_success,"
             << "modify_maxsat_time,modify_modif,modify_a,modify_b,modify_success\n";
        
        // Data
        for (const auto& r : results) {
            file << r.m << "," << r.n_objects << "," << r.k << "," << r.exp_id << ","
                 << r.user_maxsat_time << "," << r.user_modif << "," << r.user_a << "," 
                 << r.user_b << "," << r.user_success << ","
                 << r.resource_maxsat_time << "," << r.resource_modif << "," << r.resource_a << "," 
                 << r.resource_b << "," << r.resource_success << ","
                 << r.modify_maxsat_time << "," << r.modify_modif << "," << r.modify_a << "," 
                 << r.modify_b << "," << r.modify_success << "\n";
        }
        
        file.close();
        std::cout << "" << results.size() << " résultats sauvegardés dans " << output_file << std::endl;
        return true;
    }
    
    /**
     * @brief Charge des résultats existants depuis un CSV
     */
    bool load() {
        std::ifstream file(output_file);
        if (!file.is_open()) return false;
        
        std::string line;
        std::getline(file, line);  // Skip header
        
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string val;
            ExperimentResult r;
            
            std::getline(ss, val, ','); r.m = std::stoi(val);
            std::getline(ss, val, ','); r.n_objects = std::stoi(val);
            std::getline(ss, val, ','); r.k = std::stoi(val);
            std::getline(ss, val, ','); r.exp_id = std::stoi(val);
            std::getline(ss, val, ','); r.user_maxsat_time = std::stod(val);
            std::getline(ss, val, ','); r.user_modif = std::stoi(val);
            std::getline(ss, val, ','); r.user_a = std::stoi(val);
            std::getline(ss, val, ','); r.user_b = std::stoi(val);
            std::getline(ss, val, ','); r.user_success = std::stoi(val);
            std::getline(ss, val, ','); r.resource_maxsat_time = std::stod(val);
            std::getline(ss, val, ','); r.resource_modif = std::stoi(val);
            std::getline(ss, val, ','); r.resource_a = std::stoi(val);
            std::getline(ss, val, ','); r.resource_b = std::stoi(val);
            std::getline(ss, val, ','); r.resource_success = std::stoi(val);
            std::getline(ss, val, ','); r.modify_maxsat_time = std::stod(val);
            std::getline(ss, val, ','); r.modify_modif = std::stoi(val);
            std::getline(ss, val, ','); r.modify_a = std::stoi(val);
            std::getline(ss, val, ','); r.modify_b = std::stoi(val);
            std::getline(ss, val, ','); r.modify_success = std::stoi(val);
            
            results.push_back(r);
        }
        
        file.close();
        std::cout << "" << results.size() << " résultats chargés depuis " << output_file << std::endl;
        return true;
    }
    
    void clear() { results.clear(); }
    size_t size() const { return results.size(); }
};
