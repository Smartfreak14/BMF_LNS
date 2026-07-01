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
