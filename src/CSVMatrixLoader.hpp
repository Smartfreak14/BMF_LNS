/**
 * @file CSVMatrixLoader.hpp
 * @brief Loader simple pour matrices binaires CSV (sans dépendance MaLib)
 *        Parsing direct du format CSV avec support des valeurs manquantes
 */

#ifndef CSV_MATRIX_LOADER_HPP
#define CSV_MATRIX_LOADER_HPP

#include "Matrix.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <algorithm>

/**
 * @class CSVMatrixLoader
 * @brief Charge des matrices binaires depuis des fichiers CSV
 */
class CSVMatrixLoader {
public:
    /**
     * @brief Charge une matrice depuis un fichier CSV
     * @param filename Chemin vers le fichier CSV
     * @return Matrix chargée (les valeurs manquantes sont marquées comme -1)
     */
    static Matrix loadFromCSV(const std::string& filename) {
        if (!std::filesystem::exists(filename)) {
            std::cerr << "Erreur: fichier non trouvé: " << filename << std::endl;
            return Matrix(0, 0);
        }
        
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Erreur: impossible d'ouvrir " << filename << std::endl;
            return Matrix(0, 0);
        }
        
        // Première passe: compter les dimensions
        std::vector<std::vector<int>> data;
        std::string line;
        
        while (std::getline(file, line)) {
            // Ignorer les lignes vides
            if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos) {
                continue;
            }
            
            std::vector<int> row;
            std::stringstream ss(line);
            std::string cell;
            
            while (std::getline(ss, cell, ',')) {
                // Nettoyer la cellule
                cell.erase(std::remove_if(cell.begin(), cell.end(), 
                    [](char c) { return c == ' ' || c == '\t' || c == '\r'; }), cell.end());
                
                if (cell.empty() || cell == "?") {
                    row.push_back(-1);  // Valeur manquante
                } else {
                    try {
                        int val = std::stoi(cell);
                        row.push_back(val == 0 ? 0 : 1);  // Binaire
                    } catch (...) {
                        row.push_back(-1);  // Valeur non-numérique = manquante
                    }
                }
            }
            
            if (!row.empty()) {
                data.push_back(row);
            }
        }
        file.close();
        
        if (data.empty()) {
            std::cerr << "Erreur: matrice vide dans " << filename << std::endl;
            return Matrix(0, 0);
        }
        
        int m = data.size();
        int n = data[0].size();
        
        // Vérifier la cohérence des colonnes
        for (const auto& row : data) {
            if ((int)row.size() != n) {
                n = std::max(n, (int)row.size());
            }
        }
        
        // Créer la matrice avec les valeurs
        Matrix M(m, n, 0);
        
        int ones = 0, zeros = 0, missing = 0;
        
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < (int)data[i].size() && j < n; j++) {
                M(i, j) = data[i][j];
                if (data[i][j] == 1) ones++;
                else if (data[i][j] == 0) zeros++;
                else missing++;
            }
        }
        
        std::cout << "  Matrice chargée: " << m << "x" << n << std::endl;
        std::cout << "  - Ones: " << ones << " (" << (100.0 * ones / (m * n)) << "%)" << std::endl;
        std::cout << "  - Zeros: " << zeros << " (" << (100.0 * zeros / (m * n)) << "%)" << std::endl;
        if (missing > 0) {
            std::cout << "  - Missing: " << missing << " (" << (100.0 * missing / (m * n)) << "%)" << std::endl;
        }
        
        return M;
    }
    
    /**
     * @brief Liste tous les fichiers CSV dans un répertoire
     * @param directory Chemin vers le répertoire
     * @return Vecteur des chemins vers les fichiers CSV
     */
    static std::vector<std::string> listCSVFiles(const std::string& directory) {
        std::vector<std::string> files;
        
        if (!std::filesystem::exists(directory)) {
            std::cerr << "Erreur: répertoire non trouvé: " << directory << std::endl;
            return files;
        }
        
        for (const auto& entry : std::filesystem::directory_iterator(directory)) {
            if (entry.path().extension() == ".csv") {
                files.push_back(entry.path().string());
            }
        }
        
        std::sort(files.begin(), files.end());
        return files;
    }
    
    /**
     * @brief Obtient le nom de base d'un fichier (sans chemin ni extension)
     */
    static std::string getBaseName(const std::string& filepath) {
        return std::filesystem::path(filepath).stem().string();
    }
};

#endif // CSV_MATRIX_LOADER_HPP
