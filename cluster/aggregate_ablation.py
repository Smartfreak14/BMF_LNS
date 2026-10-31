#!/usr/bin/env python3
"""
Agrege les CSV d'ablation (4 configurations par dataset, multi-seed possible)
en un CSV wide-format pour le tableau d'ablation du papier.

Reconnait DEUX patterns de nom de fichier :
  1. Ancien : <dataset>_<mode>.csv          (single seed, 1 CSV par mode)
  2. Nouveau : <dataset>_<mode>_s<N>.csv    (multi-seed, N CSVs par mode)

Si multi-seed : renvoie la MEDIANE de err et time_ms sur les seeds, plus
min_err (best case) et err_std. Colonnes de sortie :

  dataset, k, m, n, regime,
  greedy_err_median, greedy_err_min, greedy_time_ms_median, greedy_nseeds,
  no_double_flip_err_median, ..., no_double_flip_nseeds,
  no_lns_err_median, ..., no_lns_nseeds,
  full_err_median, ..., full_nseeds

Usage :
    python aggregate_ablation.py results_ablation_real/ -o ablation_all.csv
"""
import argparse, csv, re, sys, statistics
from pathlib import Path

MODES = ["greedy", "no-weighting", "wls-no-double-no-lns",
         "no-double-flip", "no-lns", "full"]
MODE_SAFE = {"greedy": "greedy",
             "no-weighting": "no_weighting",
             "wls-no-double-no-lns": "wls_no_df_no_lns",
             "no-double-flip": "no_double_flip",
             "no-lns": "no_lns",
             "full": "full"}

# Regex pour capturer <dataset>_<mode>[_s<N>].csv (les modes plus longs en 1er)
PAT_MULTI = re.compile(
    r"^(?P<dataset>.+?)_(?P<mode>wls-no-double-no-lns|no-double-flip|no-weighting|greedy|no-lns|full)(?:_s(?P<seed>\d+))?$"
)


def parse_stem(stem: str):
    m = PAT_MULTI.match(stem)
    if not m:
        return None, None, None
    return m.group("dataset") + ".csv", m.group("mode"), m.group("seed")


def collect(results_dir: Path):
    # (dataset, mode) -> list of (err, time_ms, phase1_err, verif, m, n, regime, k)
    groups = {}
    for csv_path in sorted(results_dir.glob("*.csv")):
        dataset, mode, seed = parse_stem(csv_path.stem)
        if dataset is None or mode is None:
            continue
        with csv_path.open() as f:
            reader = csv.DictReader(f)
            for r in reader:
                if not r.get("dataset"):
                    continue
                try:
                    err = int(r["final_errors"])
                    tms = float(r["time_ms"])
                except (KeyError, ValueError, TypeError):
                    continue
                p1_err = r.get("phase1_errors", "")
                verif = r.get("verif_errors", "")
                m = r.get("m", "")
                n = r.get("n", "")
                regime = r.get("regime", "")
                k = r.get("k", "")
                key = (dataset, mode)
                groups.setdefault(key, {"runs": [], "k": k, "m": m, "n": n, "regime": regime})
                groups[key]["runs"].append({
                    "err": err, "time_ms": tms,
                    "phase1_err": p1_err, "verif": verif,
                })
    return groups


def summarize(runs):
    errs = [r["err"] for r in runs]
    times = [r["time_ms"] for r in runs]
    return {
        "err_median": int(statistics.median(errs)),
        "err_min": min(errs),
        "err_max": max(errs),
        "err_std": round(statistics.pstdev(errs), 2) if len(errs) > 1 else 0,
        "time_ms_median": round(statistics.median(times), 1),
        "nseeds": len(runs),
        "verif": runs[0].get("verif", ""),
        "phase1_err": runs[0].get("phase1_err", ""),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("results_dir", type=Path)
    p.add_argument("-o", "--output", type=Path, default=Path("ablation_all.csv"))
    args = p.parse_args()

    if not args.results_dir.is_dir():
        print(f"[error] {args.results_dir} n'est pas un dossier", file=sys.stderr)
        return 1

    groups = collect(args.results_dir)
    if not groups:
        print(f"[error] aucun CSV exploitable dans {args.results_dir}", file=sys.stderr)
        print(f"        pattern attendu : <dataset>_<mode>[_s<N>].csv",
              file=sys.stderr)
        print(f"        modes reconnus : {', '.join(MODES)}", file=sys.stderr)
        return 1

    # (dataset, mode) -> summary
    dataset_meta = {}  # dataset -> (k, m, n, regime)
    dataset_mode = {}  # (dataset, mode) -> summary
    for (dataset, mode), g in groups.items():
        dataset_meta.setdefault(dataset, (g["k"], g["m"], g["n"], g["regime"]))
        dataset_mode[(dataset, mode)] = summarize(g["runs"])

    header = ["dataset", "k", "m", "n", "regime"]
    # Ordre canonique (cumulative) : greedy -> no-weighting -> WLS-only -> +DF -> full
    for m in ["greedy", "no_weighting", "wls_no_df_no_lns",
              "no_double_flip", "no_lns", "full"]:
        header += [m + "_err_median", m + "_err_min", m + "_err_std",
                   m + "_time_ms_median", m + "_phase1_err", m + "_verif",
                   m + "_nseeds"]

    with args.output.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction="ignore")
        w.writeheader()
        n_full = 0
        for dataset in sorted(dataset_meta):
            k, mm, nn, regime = dataset_meta[dataset]
            row = {"dataset": dataset, "k": k, "m": mm, "n": nn, "regime": regime}
            covered_modes = 0
            for mode in MODES:
                pref = MODE_SAFE[mode]
                summ = dataset_mode.get((dataset, mode))
                if summ is None:
                    continue
                row[pref + "_err_median"] = summ["err_median"]
                row[pref + "_err_min"] = summ["err_min"]
                row[pref + "_err_std"] = summ["err_std"]
                row[pref + "_time_ms_median"] = summ["time_ms_median"]
                row[pref + "_phase1_err"] = summ["phase1_err"]
                row[pref + "_verif"] = summ["verif"]
                row[pref + "_nseeds"] = summ["nseeds"]
                covered_modes += 1
            w.writerow(row)
            if covered_modes == 5:
                n_full += 1

    total = len(dataset_meta)
    print(f"[ok] {total} datasets ecrits ({n_full} avec les 5 configs completes) dans {args.output}")
    if n_full < total:
        print(f"[warn] {total - n_full} dataset(s) ont moins de 5 configs -> voir les runs manquants",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
