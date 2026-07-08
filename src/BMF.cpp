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
#include <iomanip>


void BMF::dump_wcnf(const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out.is_open()) {
        std::cerr << "dump_wcnf: cannot open " << output_path << std::endl;
        return;
    }

    // A[i,l] : id = 1 + i*k + l
    // B[l,j] : id = 1 + m*k + l*n + j
    // p[i,l,j] : id auxiliaire consecutif a partir de m*k+k*n+1 (cellules M=1)
    auto A_var = [&](int i, int l) { return 1 + i * k + l; };
    auto B_var = [&](int l, int j) { return 1 + m * k + l * n + j; };

    // Passe 1 : compter clauses/vars/top pour le header v1.
    long long n_hard = 0, n_soft = 0, total_soft_weight = 0;
    long long n_M1_cells = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int m_ij = M(i, j);
            if (m_ij == -1) continue;
            if (m_ij == 1) {
                n_hard += 3 * (long long)k;
                n_soft += 1;
                total_soft_weight += 1;
                n_M1_cells++;
            } else {
                n_hard += (long long)k;
            }
        }
    }
    long long n_aux = n_M1_cells * (long long)k;
    long long total_vars = (long long)m * k + (long long)k * n + n_aux;
    long long top = total_soft_weight + 1;

    out << "c BMF: m=" << m << " n=" << n << " k=" << k
        << " M1_cells=" << n_M1_cells << "\n"
        << "c Encoding: M=1 -> soft (∃l: A[i,l] ∧ B[l,j]) via Tseitin p_{i,l,j}\n"
        << "c           M=0 -> hard ¬A[i,l] ∨ ¬B[l,j] pour chaque l\n"
        << "c Vars: A[i,l]=1+i*k+l, B[l,j]=1+m*k+l*n+j, aux p_{i,l,j}\n"
        << "p wcnf " << total_vars << " " << (n_hard + n_soft) << " " << top << "\n";

    // Passe 2 : ecrire les clauses.
    int next_aux = 1 + m * k + k * n;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int m_ij = M(i, j);
            if (m_ij == -1) continue;

            if (m_ij == 1) {
                std::vector<int> big_or;
                big_or.reserve(k);
                for (int l = 0; l < k; l++) {
                    int p = next_aux++;
                    int a = A_var(i, l);
                    int b = B_var(l, j);
                    out << top << ' ' << -p << ' ' <<  a << " 0\n";
                    out << top << ' ' << -p << ' ' <<  b << " 0\n";
                    out << top << ' ' << -a << ' ' << -b << ' ' << p << " 0\n";
                    big_or.push_back(p);
                }
                out << 1;
                for (int lit : big_or) out << ' ' << lit;
                out << " 0\n";
            } else {
                for (int l = 0; l < k; l++) {
                    out << top << ' ' << -A_var(i, l) << ' ' << -B_var(l, j) << " 0\n";
                }
            }
        }
    }

    std::cout << "[wcnf] " << output_path
              << " nvars=" << total_vars
              << " nclauses=" << (n_hard + n_soft)
              << " (hard=" << n_hard << " soft=" << n_soft << ")"
              << " top=" << top << std::endl;
}


int BMF::lns_step_partial(
    Matrix& A_current, Matrix& B_current,
    const std::vector<int>& rows_to_free,
    const std::vector<int>& cols_to_free
) {
    std::vector<bool> free_rows(m, false);
    std::vector<bool> free_cols(n, false);

    for (int r : rows_to_free) if (r >= 0 && r < m) free_rows[r] = true;
    for (int c : cols_to_free) if (c >= 0 && c < n) free_cols[c] = true;

    int fixed_cost = 0;
    Matrix computed = A_current.multiply(B_current);

    for (int i = 0; i < m; i++) {
        if (free_rows[i]) continue;
        for (int j = 0; j < n; j++) {
            if (free_cols[j]) continue;
            if (M(i, j) == -1) continue;
            if (computed(i, j) != M(i, j))
                fixed_cost++;
        }
    }

    SATSolver partial_solver;
    VariableManager A_partial(partial_solver), B_partial(partial_solver);

    for (int i = 0; i < m; i++)
        if (free_rows[i])
            for (int l = 0; l < k; l++)
                A_partial.get(i, l);

    for (int l = 0; l < k; l++)
        for (int j = 0; j < n; j++)
            if (free_cols[j])
                B_partial.get(l, j);

    int hard_count = 0, soft_count = 0;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            bool row_free = free_rows[i];
            bool col_free = free_cols[j];

            if (!row_free && !col_free) continue;
            if (M(i, j) == -1) continue;

            bool is_error = (computed(i, j) != M(i, j));

            if (M(i, j) == 1) {
                std::vector<int> big_or;

                for (int l = 0; l < k; l++) {
                    if (!row_free) {
                        if (A_current(i, l) == 1)
                            big_or.push_back(B_partial.get(l, j));
                    } else if (!col_free) {
                        if (B_current(l, j) == 1)
                            big_or.push_back(A_partial.get(i, l));
                    } else {
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
                    if (is_error) {
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

            } else {
                for (int l = 0; l < k; l++) {
                    std::vector<int> clause;

                    if (!row_free) {
                        if (A_current(i, l) == 1)
                            clause = {-B_partial.get(l, j)};
                    } else if (!col_free) {
                        if (B_current(l, j) == 1)
                            clause = {-A_partial.get(i, l)};
                    } else {
                        clause = {-A_partial.get(i, l), -B_partial.get(l, j)};
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

    std::vector<int> assumptions;

    for (int i = 0; i < m; i++) {
        if (free_rows[i]) continue;
        for (int l = 0; l < k; l++) {
            auto it = A_partial.vars.find({i, l});
            if (it != A_partial.vars.end())
                assumptions.push_back(A_current(i, l) == 1 ? it->second : -it->second);
        }
    }

    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            if (free_cols[j]) continue;
            auto it = B_partial.vars.find({l, j});
            if (it != B_partial.vars.end())
                assumptions.push_back(B_current(l, j) == 1 ? it->second : -it->second);
        }
    }

    if (!partial_solver.solve_with_assumptions(assumptions))
        return -1;

    long long partial_cost = partial_solver.get_cost();

    for (int i = 0; i < m; i++) {
        if (!free_rows[i]) continue;
        for (int l = 0; l < k; l++) {
            auto it = A_partial.vars.find({i, l});
            if (it != A_partial.vars.end())
                A_current(i, l) = partial_solver.getValue(it->second) ? 1 : 0;
        }
    }

    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            if (!free_cols[j]) continue;
            auto it = B_partial.vars.find({l, j});
            if (it != B_partial.vars.end())
                B_current(l, j) = partial_solver.getValue(it->second) ? 1 : 0;
        }
    }

    return fixed_cost + static_cast<int>(partial_cost);
}
