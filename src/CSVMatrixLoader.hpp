#pragma once
#include "Matrix.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <algorithm>

/**
 * @class CSVMatrixLoader
 * @brief Loads binary matrices from CSV files.
 */
class CSVMatrixLoader {
public:
    /**
     * @brief Loads a matrix from a CSV file.
     * @param filename Path to the CSV file.
     * @return Loaded matrix (missing values are marked as -1).
     */
    static Matrix loadFromCSV(const std::string& filename) {
        if (!std::filesystem::exists(filename)) {
            std::cerr << "Error: file not found: " << filename << std::endl;
            return Matrix(0, 0);
        }

        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: cannot open " << filename << std::endl;
            return Matrix(0, 0);
        }

        std::vector<std::vector<int>> data;
        std::string line;

        while (std::getline(file, line)) {
            if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos) {
                continue;
            }

            std::vector<int> row;
            std::stringstream ss(line);
            std::string cell;

            while (std::getline(ss, cell, ',')) {
                cell.erase(std::remove_if(cell.begin(), cell.end(),
                    [](char c) { return c == ' ' || c == '\t' || c == '\r'; }), cell.end());

                if (cell.empty() || cell == "?") {
                    row.push_back(-1);
                } else {
                    try {
                        int val = std::stoi(cell);
                        row.push_back(val == 0 ? 0 : 1);
                    } catch (...) {
                        row.push_back(-1);
                    }
                }
            }

            if (!row.empty()) {
                data.push_back(row);
            }
        }
        file.close();

        if (data.empty()) {
            std::cerr << "Error: empty matrix in " << filename << std::endl;
            return Matrix(0, 0);
        }

        int m = data.size();
        int n = data[0].size();

        for (const auto& row : data) {
            if ((int)row.size() != n) {
                n = std::max(n, (int)row.size());
            }
        }

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

        std::cout << "  Matrix loaded: " << m << "x" << n << std::endl;
        std::cout << "  - Ones: " << ones << " (" << (100.0 * ones / (m * n)) << "%)" << std::endl;
        std::cout << "  - Zeros: " << zeros << " (" << (100.0 * zeros / (m * n)) << "%)" << std::endl;
        if (missing > 0) {
            std::cout << "  - Missing: " << missing << " (" << (100.0 * missing / (m * n)) << "%)" << std::endl;
        }

        return M;
    }

    /**
     * @brief Lists all CSV files in a directory.
     * @param directory Path to the directory.
     * @return Vector of paths to CSV files.
     */
    static std::vector<std::string> listCSVFiles(const std::string& directory) {
        std::vector<std::string> files;

        if (!std::filesystem::exists(directory)) {
            std::cerr << "Error: directory not found: " << directory << std::endl;
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

    /** @brief Gets the base name of a file (without path or extension). */
    static std::string getBaseName(const std::string& filepath) {
        return std::filesystem::path(filepath).stem().string();
    }
};
