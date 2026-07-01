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
};
