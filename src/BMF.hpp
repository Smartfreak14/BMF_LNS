#pragma once
#include "Matrix.hpp"
#include "SATSolver.hpp"
#include <vector>

class BMF {
public:
    int k, m, n;
    Matrix M;
    SATSolver solver;
    VariableManager A, B;

    BMF(int k_factors, const Matrix& matrix)
        : k(k_factors), M(matrix), A(solver), B(solver) {
        m = M.rows;
        n = M.cols;
    }

    int lns_step_partial(
        Matrix& A_current, Matrix& B_current,
        const std::vector<int>& rows_to_free,
        const std::vector<int>& cols_to_free
    );

    // Ecrit le probleme BMF complet au format WCNF v1 (Partial MaxSAT).
    // Encodage : M=1 -> soft (exists l: A[i,l] AND B[l,j]) via Tseitin ;
    //            M=0 -> hard ¬A[i,l] ∨ ¬B[l,j] pour chaque l ;
    //            M=-1 -> ignore.
    void dump_wcnf(const std::string& output_path);
};
