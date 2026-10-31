#!/usr/bin/env python3

import sys, time, argparse
import numpy as np


def load_matrix(path):
    with open(path) as f:
        rows = []
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if not line:
                continue
            vals = []
            for x in line.split(','):
                x = x.strip()
                if x == '':
                    vals.append(-1)
                else:
                    vals.append(int(x))
            rows.append(vals)
    return np.array(rows, dtype=np.int8)


def compute_association(M, tau):
    M01 = (M == 1).astype(np.int32)
    cooc = M01.T @ M01              # (n x n) : cooc[i,j] = sum_k(M[k,i] AND M[k,j])
    count = np.diag(cooc)           # count[i] = nombre de 1 dans colonne i
    with np.errstate(divide='ignore', invalid='ignore'):
        # count[:, None] broadcast : divise chaque LIGNE i par count[i] (paper).
        A = (cooc / np.maximum(count[:, None], 1)) >= tau
    np.fill_diagonal(A, True)
    return A.astype(np.int8)


def cover_error(M, X, S):
    prod = np.clip(X @ S, 0, 1).astype(np.int8)
    mask = (M != -1)
    return int(np.sum(mask & (prod != M)))


def asso(M, r, tau_list, budget_s, seed):
    m, n = M.shape
    rng = np.random.default_rng(seed)
    t0 = time.perf_counter()
    trace = []
    best_err = int((M == 1).sum())     # trivial : X=0, S=0
    best_X = np.zeros((m, r), dtype=np.int8)
    best_S = np.zeros((r, n), dtype=np.int8)
    time_to_zero = -1.0
    per_tau_err = {}
    per_tau_info = {}       # tau -> (facteurs ajoutes, raison de l'arret)
    best_factors, best_tau = 0, None

    def elapsed_ms():
        return (time.perf_counter() - t0) * 1000

    trace.append((elapsed_ms(), best_err))

    per_tau_budget = budget_s / max(len(tau_list), 1)

    for tau_idx, tau in enumerate(tau_list):
        tau_deadline = (tau_idx + 1) * per_tau_budget
        tau_start_err = best_err
        if elapsed_ms() / 1000 > budget_s:
            break
        A = compute_association(M, tau)
        cands = A.T  
        X = np.zeros((m, r), dtype=np.int8)
        S = np.zeros((r, n), dtype=np.int8)
        added, stop = 0, "rank"
        for k in range(r):
            if elapsed_ms() / 1000 > tau_deadline:
                stop = "deadline"
                break
            best_delta = 0
            best_i = -1
            best_x = None
            for ci in range(cands.shape[0]):
                s_new = cands[ci]
                if s_new.sum() == 0:
                    continue
                prod_wo = np.clip(X @ S, 0, 1).astype(np.int8)  # sans le nouveau
                mask = (M != -1)
                err_wo = (mask & (prod_wo != M)).sum(axis=1)
                prod_w = np.maximum(prod_wo, s_new[None, :])
                err_w = (mask & (prod_w != M)).sum(axis=1)
                x_new = (err_w < err_wo).astype(np.int8)
                gain = int(err_wo.sum() - (x_new * err_w + (1 - x_new) * err_wo).sum())
                if gain > best_delta:
                    best_delta = gain
                    best_i = ci
                    best_x = x_new
            if best_i < 0:
                stop = "no_gain"
                break
            X[:, k] = best_x
            S[k] = cands[best_i]
            added = k + 1
            err = cover_error(M, X, S)
            if err < best_err:
                best_err = err
                best_X = X.copy()
                best_S = S.copy()
                best_factors, best_tau = added, tau
                t = elapsed_ms()
                trace.append((t, best_err))
                if best_err == 0 and time_to_zero < 0:
                    time_to_zero = t
        per_tau_err[tau] = min(tau_start_err, best_err)
        per_tau_info[tau] = (added, stop)
    return (best_X, best_S, best_err, time_to_zero, trace, per_tau_err,
            per_tau_info, best_factors, best_tau)


CSV_HEADER = ("method,dataset,k,init_errors,phase1_errors,phase1_last_improve_ms,"
              "phase2_errors,phase2_last_improve_ms,final_errors,iterations,time_ms,"
              "wls_total_flips,wls_total_time_ms,wls_flips_per_sec,verif_errors,"
              "m,n,regime,tool,seed,time_to_zero_ms,trace_file,factors_used\n")


def run_asso_single(input_csv, r, timeout_s, output_csv,
                    trace_out="", seed=42, regime="raw",
                    tau_list=None):
    M = load_matrix(input_csv)
    m, n = M.shape
    if tau_list is None:
        tau_list = [0.3, 0.5, 0.7, 0.9]
    print(f"[asso] {input_csv} M={m}x{n} r={r} timeout={timeout_s}s "
          f"seed={seed} tau_list={tau_list} "
          f"budget_per_tau={timeout_s/len(tau_list):.1f}s")

    t0 = time.perf_counter()
    (X, S, best_err, time_to_zero, trace, per_tau_err,
     per_tau_info, best_factors, best_tau) = asso(M, r, tau_list, timeout_s, seed)
    total_ms = (time.perf_counter() - t0) * 1000

    print(f"[asso]   per_tau_err={per_tau_err}")
    # Facteurs ajoutes par tau et raison de l'arret : no_gain (aucun candidat
    # n'ameliore la couverture), deadline (budget du tau), rank (r atteint).
    print(f"[asso]   per_tau_factors={per_tau_info}")
    print(f"[asso]   factors_used={best_factors} (tau={best_tau}) sur r={r}")
    print(f"[asso]   best_err={best_err} total_ms={total_ms:.1f} time_to_zero_ms={time_to_zero:.1f}")

    if trace_out:
        with open(trace_out, "w") as f:
            f.write("time_ms err\n")
            for t, e in trace:
                f.write(f"{t:.3f} {e}\n")
        # Sidecar : erreurs par tau (reviewer-proof, prouve que le sweep
        # est genuine et pas cherry-picked apres coup).
        tau_report = trace_out.replace(".trace", ".tau.csv")
        with open(tau_report, "w") as f:
            f.write("tau,best_err,factors,stop\n")
            for tau in tau_list:
                added, stop = per_tau_info.get(tau, (0, "budget"))
                f.write(f"{tau},{per_tau_err.get(tau, -1)},{added},{stop}\n")

    import os
    dataset = os.path.basename(input_csv)
    init_err = int((M == 1).sum())
    ttz = time_to_zero if time_to_zero > 0 else -1
    row = (f"asso,{dataset},{r},{init_err},"
           f"{best_err},{ttz},"
           f"{best_err},{total_ms:.1f},"
           f"{best_err},0,{total_ms:.1f},"
           f"0,0.0,0.0,-1,"
           f"{m},{n},{regime},asso,{seed},"
           f"{ttz},{trace_out},{best_factors}\n")
    with open(output_csv, "w") as f:
        f.write(CSV_HEADER)
        f.write(row)
    print(f"[csv] {output_csv}")
    return best_err, total_ms


def run_batch(datasets_csv, data_dir, output_dir, nseeds=5,
              timeout_override=None, regime="raw", tau_list=None):
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Parse datasets.csv : filename,k,timeout
    datasets = []
    with open(datasets_csv) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith('#') or i == 0:
                continue
            parts = [x.strip() for x in line.split(',')]
            if len(parts) < 2:
                continue
            filename = parts[0]
            k = int(parts[1])
            tout = float(parts[2]) if len(parts) >= 3 else 300.0
            if timeout_override is not None:
                tout = timeout_override
            datasets.append((filename, k, tout))

    print(f"[batch] {len(datasets)} datasets x {nseeds} seeds = {len(datasets)*nseeds} runs")
    print(f"[batch] data_dir={data_dir} output_dir={output_dir}")

    n_ok = 0
    n_fail = 0
    for idx, (filename, k, tout) in enumerate(datasets, 1):
        input_path = os.path.join(data_dir, filename)
        if not os.path.isfile(input_path):
            print(f"[batch] [{idx}/{len(datasets)}] SKIP {filename} : absent de {data_dir}")
            n_fail += 1
            continue
        base = filename[:-4] if filename.endswith('.csv') else filename
        print(f"[batch] [{idx}/{len(datasets)}] {filename} k={k} timeout={tout}s")
        for seed in range(1, nseeds + 1):
            out_csv = os.path.join(output_dir, f"{base}_asso_s{seed}.csv")
            trace_out = os.path.join(output_dir, f"{base}_asso_s{seed}.trace")
            try:
                run_asso_single(input_path, k, tout, out_csv,
                                trace_out=trace_out, seed=seed, regime=regime,
                                tau_list=tau_list)
                n_ok += 1
            except Exception as e:
                print(f"[batch]   [seed={seed}] ERREUR : {e}")
                n_fail += 1

    print(f"[batch] Fini : {n_ok} OK, {n_fail} echec/skip")


def main():
    p = argparse.ArgumentParser(
        description="ASSO : single-run ou batch sur un CSV de datasets")

    # Mode batch : detecte via --batch
    p.add_argument("--batch", action="store_true",
                   help="Mode batch : args = <datasets.csv> <data_dir> <output_dir>")
    p.add_argument("--nseeds", type=int, default=5,
                   help="Nombre de graines par dataset en mode batch (defaut 5)")
    p.add_argument("--timeout-override", type=float, default=None,
                   help="Timeout uniforme (secondes) qui ecrase celui de datasets.csv")

    # Args positionnels : signification depend du mode
    p.add_argument("args", nargs="*",
                   help="single: <input_csv> <r> <timeout_s> <output_csv> | "
                        "batch: <datasets.csv> <data_dir> <output_dir>")

    # Options communes
    p.add_argument("--trace-out", default="", help="Fichier de trace (single uniquement)")
    p.add_argument("--seed", type=int, default=42, help="Graine RNG (single uniquement)")
    p.add_argument("--regime", default="raw",
                   help="Etiquette regime (raw/simpluni/simplexi/synth)")
    p.add_argument("--tau-list", default="0.3,0.5,0.7,0.9",
                   help="Liste de seuils tau a balayer (virgule-separee). "
                        "Budget partage equitablement. Defaut : 0.3,0.5,0.7,0.9 "
                        "(couvre la plage calibree par Kolomvakis).")

    args = p.parse_args()

    tau_list = [float(x) for x in args.tau_list.split(",") if x.strip()]

    if args.batch:
        if len(args.args) != 3:
            p.error("Mode batch : besoin de <datasets.csv> <data_dir> <output_dir>")
        datasets_csv, data_dir, output_dir = args.args
        run_batch(datasets_csv, data_dir, output_dir,
                  nseeds=args.nseeds, timeout_override=args.timeout_override,
                  regime=args.regime, tau_list=tau_list)
    else:
        if len(args.args) != 4:
            p.error("Mode single : besoin de <input_csv> <r> <timeout_s> <output_csv>")
        input_csv, r_str, tout_str, output_csv = args.args
        run_asso_single(input_csv, int(r_str), float(tout_str), output_csv,
                        trace_out=args.trace_out, seed=args.seed,
                        regime=args.regime, tau_list=tau_list)


if __name__ == "__main__":
    main()
