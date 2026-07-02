#include "BMFLocalSearch.hpp"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <climits>
#include <limits>

BMFLocalSearch::BMFLocalSearch(int k_factors, const Matrix& matrix, unsigned seed)
    : m(matrix.rows), n(matrix.cols), k(k_factors), M(matrix),
      A(m, k, 0), B(k, n, 0),
      rng(seed) {}


void BMFLocalSearch::initialize_greedy() {
    std::vector<std::vector<bool>> uncovered(m, std::vector<bool>(n, false));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            if (M(i, j) == 1)
                uncovered[i][j] = true;

    for (int i = 0; i < m; i++)
        for (int l = 0; l < k; l++)
            A(i, l) = 0;
    for (int l = 0; l < k; l++)
        for (int j = 0; j < n; j++)
            B(l, j) = 0;

    for (int l = 0; l < k; l++) {
        int best_col = -1;
        int best_count = 0;
        for (int j = 0; j < n; j++) {
            int cnt = 0;
            for (int i = 0; i < m; i++)
                if (uncovered[i][j]) cnt++;
            if (cnt > best_count) { best_count = cnt; best_col = j; }
        }
        if (best_col == -1 || best_count == 0) break;

        B(l, best_col) = 1;
        for (int i = 0; i < m; i++) {
            if (M(i, best_col) == 1) {
                A(i, l) = 1;
                for (int j = 0; j < n; j++)
                    if (B(l, j) == 1 && uncovered[i][j])
                        uncovered[i][j] = false;
            }
        }

        for (int j = 0; j < n; j++) {
            if (B(l, j) == 0) {
                int gain = 0;
                for (int i = 0; i < m; i++) {
                    if (A(i, l) == 1) {
                        if (uncovered[i][j]) gain++;
                        else if (M(i, j) == 0) gain--;
                    }
                }
                if (gain > 0) {
                    B(l, j) = 1;
                    for (int i = 0; i < m; i++)
                        if (A(i, l) == 1 && uncovered[i][j])
                            uncovered[i][j] = false;
                }
            }
        }
    }
}


// Rebuild toutes les structures depuis (A, B, M). Appeler apres modification
// externe de A/B.
void BMFLocalSearch::compute_all_counts() {
    nbcover_flat.assign(m * n, 0);

    rows_A.assign(k, {});
    cols_B.assign(k, {});
    cols_A.assign(m, {});
    rows_B.assign(n, {});
    zero_cover_in_row.assign(m, {});
    zero_cover_in_col.assign(n, {});

    for (int l = 0; l < k; l++) {
        rows_A[l].resize_capacity(m);
        cols_B[l].resize_capacity(n);
    }
    for (int i = 0; i < m; i++) {
        cols_A[i].resize_capacity(k);
        zero_cover_in_row[i].resize_capacity(n);
    }
    for (int j = 0; j < n; j++) {
        rows_B[j].resize_capacity(k);
        zero_cover_in_col[j].resize_capacity(m);
    }

    score_A_flat.assign(m * k, 0.0);
    score_B_flat.assign(k * n, 0.0);
    weight_flat.assign(m * n, 1.0);

    pos_score_A.clear();
    pos_score_A.resize_capacity(m * k);
    pos_score_B.clear();
    pos_score_B.resize_capacity(k * n);
    unsat_hard_cells.clear();
    unsat_hard_cells.resize_capacity(m * n);
    unsat_soft_cells.clear();
    unsat_soft_cells.resize_capacity(m * n);

    for (int i = 0; i < m; i++) {
        for (int l = 0; l < k; l++) {
            if (A(i, l) == 1) {
                rows_A[l].insert(i);
                cols_A[i].insert(l);
            }
        }
    }
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            if (B(l, j) == 1) {
                cols_B[l].insert(j);
                rows_B[j].insert(l);
            }
        }
    }

    for (int i = 0; i < m; i++) {
        for (int l : cols_A[i].as_vector()) {
            for (int j : cols_B[l].as_vector()) {
                nbcover_flat[i * n + j]++;
            }
        }
    }

    for (int i = 0; i < m; i++) {
        int base = i * n;
        for (int j = 0; j < n; j++) {
            int nb = nbcover_flat[base + j];
            if (nb == 0) {
                zero_cover_in_row[i].insert(j);
                zero_cover_in_col[j].insert(i);
            }
            int m_ij = M(i, j);
            if (m_ij == -1) continue;
            int e = err_val(m_ij, nb);
            if (e) {
                int idx = base + j;
                if (m_ij == 1) unsat_soft_cells.insert(idx);
                else           unsat_hard_cells.insert(idx);
            }
        }
    }

    for (int i = 0; i < m; i++) {
        for (int l = 0; l < k; l++) {
            double s = 0.0;
            int sign = (A(i, l) == 1) ? -1 : +1;
            for (int j : cols_B[l].as_vector()) {
                int m_ij = M(i, j);
                if (m_ij == -1) continue;
                int nb = nbcover_flat[i * n + j];
                s += score_contrib(m_ij, nb, sign, weight_flat[i * n + j]);
            }
            score_A_flat[i * k + l] = s;
            if (s > 1e-9) pos_score_A.insert(i * k + l);
        }
    }
    for (int l = 0; l < k; l++) {
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            int sign = (B(l, j) == 1) ? -1 : +1;
            for (int i : rows_A[l].as_vector()) {
                int m_ij = M(i, j);
                if (m_ij == -1) continue;
                int nb = nbcover_flat[i * n + j];
                s += score_contrib(m_ij, nb, sign, weight_flat[i * n + j]);
            }
            score_B_flat[l * n + j] = s;
            if (s > 1e-9) pos_score_B.insert(l * n + j);
        }
    }
}


int BMFLocalSearch::count_errors() const {
    int e = 0;
    for (int i = 0; i < m; i++) {
        int base = i * n;
        for (int j = 0; j < n; j++) {
            int m_ij = M(i, j);
            if (m_ij == -1) continue;
            e += err_val(m_ij, nbcover_flat[base + j]);
        }
    }
    return e;
}


void BMFLocalSearch::flip_A(int i, int l) {
    bool old_a = (A(i, l) == 1);
    int dir = old_a ? -1 : +1;

    
    {
        // cols_B[l]/zero_cover_in_row[i] ne sont pas modifies dans cette boucle:
        // reference (pas de copie).
        const std::vector<int>& snap = zero_cover_in_row[i].as_vector();
        for (int j : snap) {
            if (B(l, j) == 1) continue;
            int m_ij = M(i, j);
            if (m_ij == -1) continue;
            double w = weight_flat[i * n + j];
            double contrib = score_contrib(m_ij, 0, +1, w);
            if (contrib == 0.0) continue;
            double delta = old_a ? -contrib : +contrib;
            score_B_flat[l * n + j] += delta;
            update_pos_score_B(l, j);
        }
    }

    // cols_B[l] n'est pas modifie par flip_A: reference (pas de copie).
    const std::vector<int>& snap_cols_B_l = cols_B[l].as_vector();

    for (int j : snap_cols_B_l) {
        int idx_ij = i * n + j;
        int old_nb = nbcover_flat[idx_ij];
        int new_nb = old_nb + dir;
        int m_ij = M(i, j);
        double w = weight_flat[idx_ij];

        const auto& rB = rows_B[j].as_vector();
        for (int l2 : rB) {
            if (l2 == l) continue;
            int sign2 = (A(i, l2) == 1) ? -1 : +1;
            double oc = score_contrib(m_ij, old_nb, sign2, w);
            double nc = score_contrib(m_ij, new_nb, sign2, w);
            double d = nc - oc;
            if (d != 0.0) {
                score_A_flat[i * k + l2] += d;
                update_pos_score_A(i, l2);
            }
        }

        const auto& cA = cols_A[i].as_vector();
        for (int l2 : cA) {
            if (l2 == l) continue;
            int sign2 = (B(l2, j) == 1) ? -1 : +1;
            double oc = score_contrib(m_ij, old_nb, sign2, w);
            double nc = score_contrib(m_ij, new_nb, sign2, w);
            double d = nc - oc;
            if (d != 0.0) {
                score_B_flat[l2 * n + j] += d;
                update_pos_score_B(l2, j);
            }
        }


        double oc_lj = old_a    ? score_contrib(m_ij, old_nb, -1, w) : 0.0;
        double nc_lj = (!old_a) ? score_contrib(m_ij, new_nb, -1, w) : 0.0;
        double d_lj = nc_lj - oc_lj;
        if (d_lj != 0.0) {
            score_B_flat[l * n + j] += d_lj;
            update_pos_score_B(l, j);
        }

        nbcover_flat[idx_ij] = new_nb;
        if (old_nb > 0 && new_nb == 0) {
            zero_cover_in_row[i].insert(j);
            zero_cover_in_col[j].insert(i);
        } else if (old_nb == 0 && new_nb > 0) {
            zero_cover_in_row[i].erase(j);
            zero_cover_in_col[j].erase(i);
        }

        if (m_ij != -1) {
            int oe = err_val(m_ij, old_nb);
            int ne = err_val(m_ij, new_nb);
            if (oe != ne) {
                if (m_ij == 1) {
                    if (ne) unsat_soft_cells.insert(idx_ij);
                    else    unsat_soft_cells.erase(idx_ij);
                } else {
                    if (ne) unsat_hard_cells.insert(idx_ij);
                    else    unsat_hard_cells.erase(idx_ij);
                }
            }
        }
    }

    A(i, l) = old_a ? 0 : 1;
    if (old_a) {
        rows_A[l].erase(i);
        cols_A[i].erase(l);
    } else {
        rows_A[l].insert(i);
        cols_A[i].insert(l);
    }

    // score_A(i, l) : direction inversee, recalcul complet sur cols_B[l].
    double s = 0.0;
    int sign = old_a ? +1 : -1;
    for (int j : snap_cols_B_l) {
        int m_ij = M(i, j);
        if (m_ij == -1) continue;
        s += score_contrib(m_ij, nbcover_flat[i * n + j], sign, weight_flat[i * n + j]);
    }
    score_A_flat[i * k + l] = s;
    update_pos_score_A(i, l);
}


void BMFLocalSearch::flip_B(int l, int j) {
    bool old_b = (B(l, j) == 1);
    int dir = old_b ? -1 : +1;

    {
        // zero_cover_in_col[j] n'est pas modifie dans cette boucle: reference.
        const std::vector<int>& snap = zero_cover_in_col[j].as_vector();
        for (int i : snap) {
            if (A(i, l) == 1) continue;
            int m_ij = M(i, j);
            if (m_ij == -1) continue;
            double w = weight_flat[i * n + j];
            double contrib = score_contrib(m_ij, 0, +1, w);
            if (contrib == 0.0) continue;
            double delta = old_b ? -contrib : +contrib;
            score_A_flat[i * k + l] += delta;
            update_pos_score_A(i, l);
        }
    }

    // rows_A[l] n'est pas modifie par flip_B: reference (pas de copie).
    const std::vector<int>& snap_rows_A_l = rows_A[l].as_vector();

    for (int i : snap_rows_A_l) {
        int idx_ij = i * n + j;
        int old_nb = nbcover_flat[idx_ij];
        int new_nb = old_nb + dir;
        int m_ij = M(i, j);
        double w = weight_flat[idx_ij];

        const auto& cA = cols_A[i].as_vector();
        for (int l2 : cA) {
            if (l2 == l) continue;
            int sign2 = (B(l2, j) == 1) ? -1 : +1;
            double oc = score_contrib(m_ij, old_nb, sign2, w);
            double nc = score_contrib(m_ij, new_nb, sign2, w);
            double d = nc - oc;
            if (d != 0.0) {
                score_B_flat[l2 * n + j] += d;
                update_pos_score_B(l2, j);
            }
        }

        const auto& rB = rows_B[j].as_vector();
        for (int l2 : rB) {
            if (l2 == l) continue;
            int sign2 = (A(i, l2) == 1) ? -1 : +1;
            double oc = score_contrib(m_ij, old_nb, sign2, w);
            double nc = score_contrib(m_ij, new_nb, sign2, w);
            double d = nc - oc;
            if (d != 0.0) {
                score_A_flat[i * k + l2] += d;
                update_pos_score_A(i, l2);
            }
        }

        double oc_il = old_b    ? score_contrib(m_ij, old_nb, -1, w) : 0.0;
        double nc_il = (!old_b) ? score_contrib(m_ij, new_nb, -1, w) : 0.0;
        double d_il = nc_il - oc_il;
        if (d_il != 0.0) {
            score_A_flat[i * k + l] += d_il;
            update_pos_score_A(i, l);
        }

        nbcover_flat[idx_ij] = new_nb;
        if (old_nb > 0 && new_nb == 0) {
            zero_cover_in_row[i].insert(j);
            zero_cover_in_col[j].insert(i);
        } else if (old_nb == 0 && new_nb > 0) {
            zero_cover_in_row[i].erase(j);
            zero_cover_in_col[j].erase(i);
        }

        if (m_ij != -1) {
            int oe = err_val(m_ij, old_nb);
            int ne = err_val(m_ij, new_nb);
            if (oe != ne) {
                if (m_ij == 1) {
                    if (ne) unsat_soft_cells.insert(idx_ij);
                    else    unsat_soft_cells.erase(idx_ij);
                } else {
                    if (ne) unsat_hard_cells.insert(idx_ij);
                    else    unsat_hard_cells.erase(idx_ij);
                }
            }
        }
    }

    B(l, j) = old_b ? 0 : 1;
    if (old_b) {
        cols_B[l].erase(j);
        rows_B[j].erase(l);
    } else {
        cols_B[l].insert(j);
        rows_B[j].insert(l);
    }

    double s = 0.0;
    int sign = old_b ? +1 : -1;
    for (int i : snap_rows_A_l) {
        int m_ij = M(i, j);
        if (m_ij == -1) continue;
        s += score_contrib(m_ij, nbcover_flat[i * n + j], sign, weight_flat[i * n + j]);
    }
    score_B_flat[l * n + j] = s;
    update_pos_score_B(l, j);
}


void BMFLocalSearch::apply_weight_change(int i, int j, double dw) {
    int idx_ij = i * n + j;
    weight_flat[idx_ij] += dw;
    int m_ij = M(i, j);
    if (m_ij == -1) return;
    int nb = nbcover_flat[idx_ij];

    for (int l : rows_B[j].as_vector()) {
        int sign = (A(i, l) == 1) ? -1 : +1;
        int de = err_val(m_ij, nb + sign) - err_val(m_ij, nb);
        if (de != 0) {
            score_A_flat[i * k + l] += -dw * (double)de;
            update_pos_score_A(i, l);
        }
    }
    for (int l : cols_A[i].as_vector()) {
        int sign = (B(l, j) == 1) ? -1 : +1;
        int de = err_val(m_ij, nb + sign) - err_val(m_ij, nb);
        if (de != 0) {
            score_B_flat[l * n + j] += -dw * (double)de;
            update_pos_score_B(l, j);
        }
    }
}


// BMS K=10 : tire K candidats aleatoires dans pos_score, garde le meilleur.
bool BMFLocalSearch::pick_good_flip(char& type, int& a, int& b) {
    int sa = (int)pos_score_A.size();
    int sb = (int)pos_score_B.size();
    int total = sa + sb;
    if (total == 0) return false;

    const int K_BMS = 10;
    double best = -std::numeric_limits<double>::infinity();
    char bt = ' ';
    int ba = -1, bb = -1;

    std::uniform_int_distribution<int> pick_pool(0, total - 1);
    for (int t = 0; t < K_BMS; t++) {
        int r = pick_pool(rng);
        if (r < sa) {
            int idx = pos_score_A.get_random(rng);
            double s = score_A_flat[idx];
            if (s > best) { best = s; bt = 'A'; ba = idx / k; bb = idx % k; }
        } else {
            int idx = pos_score_B.get_random(rng);
            double s = score_B_flat[idx];
            if (s > best) { best = s; bt = 'B'; ba = idx / n; bb = idx % n; }
        }
    }
    type = bt; a = ba; b = bb;
    return true;
}


// Cherche (i, l, j) : M=1, nbcover=0, A[i,l]=0, B[l,j]=0, score_A~0, score_B~0.
bool BMFLocalSearch::try_double_flip() {
    if (unsat_soft_cells.empty()) return false;
    const double EPS = 1e-9;

    std::vector<int> cells = unsat_soft_cells.as_vector();
    std::shuffle(cells.begin(), cells.end(), rng);

    const int MAX_CELLS = 64;
    const int MAX_K_PER_CELL = 64;
    int n_cells = std::min((int)cells.size(), MAX_CELLS);

    for (int c = 0; c < n_cells; c++) {
        int cell_idx = cells[c];
        int i = cell_idx / n;
        int j = cell_idx % n;

        std::vector<int> ls(k);
        std::iota(ls.begin(), ls.end(), 0);
        if (k > MAX_K_PER_CELL) std::shuffle(ls.begin(), ls.end(), rng);
        int n_ls = std::min(k, MAX_K_PER_CELL);

        for (int t = 0; t < n_ls; t++) {
            int l = ls[t];
            if (A(i, l) == 1) continue;
            if (B(l, j) == 1) continue;
            if (std::abs(score_A_flat[i * k + l]) > EPS) continue;
            if (std::abs(score_B_flat[l * n + j]) > EPS) continue;
            flip_A(i, l);
            flip_B(l, j);
            return true;
        }
    }
    return false;
}


// Avec proba sp : smooth (decremente satisfaites w>1). Sinon : penalise
void BMFLocalSearch::update_weights(double sp, double h_inc, double s_inc, bool& did_smooth) {
    std::uniform_real_distribution<double> u(0.0, 1.0);
    double r = u(rng);

    if (r < sp) {
        did_smooth = true;
        for (int i = 0; i < m; i++) {
            int base = i * n;
            for (int j = 0; j < n; j++) {
                int m_ij = M(i, j);
                if (m_ij == -1) continue;
                if (err_val(m_ij, nbcover_flat[base + j])) continue;
                double w = weight_flat[base + j];
                if (w <= 1.0) continue;
                double dec = (m_ij == 0) ? h_inc : s_inc;
                double new_w = std::max(1.0, w - dec);
                if (new_w != w) apply_weight_change(i, j, new_w - w);
            }
        }
        return;
    }

    did_smooth = false;
    // apply_weight_change ne change pas nbcover, donc les ensembles unsat_*
    // restent inchanges pendant l'iteration: reference (pas de copie).
    if (!unsat_hard_cells.empty()) {
        const std::vector<int>& snap = unsat_hard_cells.as_vector();
        for (int idx : snap) {
            int i = idx / n, j = idx % n;
            apply_weight_change(i, j, h_inc);
        }
    } else {
        const std::vector<int>& snap = unsat_soft_cells.as_vector();
        for (int idx : snap) {
            int i = idx / n, j = idx % n;
            apply_weight_change(i, j, s_inc);
        }
    }
}


void BMFLocalSearch::force_flip_from_unsat() {
    bool use_hard = !unsat_hard_cells.empty();
    SparseIntSet& bag = use_hard ? unsat_hard_cells : unsat_soft_cells;
    if (bag.empty()) return;
    int cell_idx = bag.get_random(rng);
    int i = cell_idx / n;
    int j = cell_idx % n;

    char best_type = ' ';
    int best_idx = -1;
    double best_score = -std::numeric_limits<double>::infinity();
    for (int l = 0; l < k; l++) {
        double sa = score_A_flat[i * k + l];
        if (sa > best_score) { best_score = sa; best_type = 'A'; best_idx = l; }
        double sb = score_B_flat[l * n + j];
        if (sb > best_score) { best_score = sb; best_type = 'B'; best_idx = l; }
    }
    if (best_type == 'A') flip_A(i, best_idx);
    else if (best_type == 'B') flip_B(best_idx, j);
}


LocalSearchResult BMFLocalSearch::solve_weighted(int max_iterations,
                                                 bool verbose,
                                                 std::atomic<bool>* stop_flag,
                                                 double timeout_ms,
                                                 double smooth_prob_override,
                                                 double h_inc,
                                                 double s_inc) {
    LocalSearchResult result;
    result.method_name = "WEIGHTED";
    auto start_time = std::chrono::high_resolution_clock::now();

    if (h_inc < 1.0) h_inc = 1.0;
    if (s_inc < 1.0) s_inc = 1.0;

    double sp = (smooth_prob_override > 0.0) ? smooth_prob_override
                                             : effective_smooth_prob();

    auto elapsed = [&]() {
        return std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - start_time).count();
    };
    auto interrupted = [&]() {
        if (stop_flag && stop_flag->load()) return true;
        if (timeout_ms > 0 && elapsed() >= timeout_ms) return true;
        return false;
    };

    compute_all_counts();

    int errors = (int)(unsat_hard_cells.size() + unsat_soft_cells.size());
    result.initial_errors = errors;
    result.error_history.push_back(errors);

    if (verbose) {
        long long nvars = (long long)m * k + (long long)k * n;
        std::cout << "  [WEIGHTED] init err=" << errors
                  << " | n_vars=" << nvars
                  << " | sp=" << std::scientific << std::setprecision(2) << sp
                  << " | h_inc=" << std::fixed << std::setprecision(1) << h_inc
                  << " | s_inc=" << std::fixed << std::setprecision(1) << s_inc
                  << std::endl;
    }

    int best_errors = errors;
    Matrix A_best = A, B_best = B;
    int iter = 0;
    long long total_flips = 0;
    int total_smooth = 0, total_penalty = 0, total_double = 0;

    while (iter < max_iterations && errors > 0 && !interrupted()) {
        char type;
        int a, b;
        if (pick_good_flip(type, a, b)) {
            if (type == 'A') flip_A(a, b);
            else             flip_B(a, b);
            total_flips++;
        } else {
            if (try_double_flip()) {
                total_double++;
                total_flips += 2;
            } else {
                bool did_smooth = false;
                update_weights(sp, h_inc, s_inc, did_smooth);
                if (did_smooth) total_smooth++; else total_penalty++;
                if (interrupted()) break;
                force_flip_from_unsat();
                total_flips++;
            }
        }

        errors = (int)(unsat_hard_cells.size() + unsat_soft_cells.size());
        result.error_history.push_back(errors);

        if (errors < best_errors) {
            best_errors = errors;
            A_best = A; B_best = B;
            result.last_improvement_ms = elapsed();
        }
        iter++;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    result.total_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    result.final_errors = best_errors;
    result.iterations = (int)total_flips;
    result.A_solution = std::move(A_best);
    result.B_solution = std::move(B_best);
    result.success = (best_errors == 0);
    result.stop_reason = (best_errors == 0) ? "zero_errors"
                       : interrupted() ? "interrupted"
                       : "max_iterations";

    if (verbose) {
        double sec = result.total_time / 1000.0;
        double fps = (sec > 0.0) ? (double)total_flips / sec : 0.0;
        std::cout << "  [WEIGHTED] Final: " << best_errors << " err in "
                  << total_flips << " flips (" << iter << " iter), "
                  << total_penalty << " penalty, " << total_smooth << " smooth, "
                  << total_double << " double-flips, "
                  << std::fixed << std::setprecision(1) << result.total_time << " ms"
                  << " (" << std::fixed << std::setprecision(0) << fps << " flips/sec)"
                  << " (" << result.stop_reason << ")" << std::endl;
    }

    return result;
}
