#pragma once
#include <vector>
#include <iostream>

// Classe de matrice binaire pour RBAC
class Matrix {
public:
    std::vector<std::vector<int>> data;
    int rows, cols;
    
    Matrix(int r, int c, int val = 0) : rows(r), cols(c) {
        data.resize(r, std::vector<int>(c, val));
    }
    
    Matrix(const std::vector<std::vector<int>>& input) {
        data = input;
        rows = data.size();
        cols = rows > 0 ? data[0].size() : 0;
    }
    
    int& operator()(int i, int j) { return data[i][j]; }
    int operator()(int i, int j) const { return data[i][j]; }
    
    // Mulutiplication booleene de matrices
    Matrix multiply(const Matrix& B) const {
        Matrix result(rows, B.cols, 0);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < B.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    if (data[i][k] == 1 && B.data[k][j] == 1) {
                        result.data[i][j] = 1;
                        break;
                    }
                }
            }
        }
        return result;
    }
    
    void print() const {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                std::cout << data[i][j] << " ";
            }
            std::cout << "\n";
        }
    }
    
    void add_row(const std::vector<int>& row) {
        data.push_back(row);
        rows++;
    }
    
    void add_cols(int n_cols, int default_val = 0) {
        for (auto& row : data) {
            for (int j = 0; j < n_cols; j++) {
                row.push_back(default_val);
            }
        }
        cols += n_cols;
    }
    
    std::vector<int>& row(int i) { return data[i]; }
    const std::vector<int>& row(int i) const { return data[i]; }
};
