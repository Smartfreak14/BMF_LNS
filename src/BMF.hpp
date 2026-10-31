#pragma once
#include "Matrix.hpp"
#include "SATSolver.hpp"
#include <atomic>
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

    // lns_step_partial dans un processus fils (EvalMaxSAT2022 ne s'interrompt pas), abandonne
    // apres limit_ms ou quand *stop devient vrai : retourne alors -2, A et B inchanges.
    int lns_step_partial_limited(
        Matrix& A_current, Matrix& B_current,
        const std::vector<int>& rows_to_free,
        const std::vector<int>& cols_to_free,
        double limit_ms, const std::atomic<bool>* stop
    );

    // Ecrit le probleme BMF complet au format WCNF v1 (Partial MaxSAT).
    // Encodage : M=1 -> soft (exists l: A[i,l] AND B[l,j]) via Tseitin ;
    //            M=0 -> hard ¬A[i,l] ∨ ¬B[l,j] pour chaque l ;
    //            M=-1 -> ignore.
    // max_clauses > 0 : si depasse, n'ecrit que l'en-tete et retourne false.
    bool dump_wcnf(const std::string& output_path, long long max_clauses = 0);

    // Ecrit l'affectation (A, B, auxiliaires p) dans la numerotation de dump_wcnf.
    void dump_assignment(const Matrix& A_val, const Matrix& B_val, const std::string& output_path);
};
