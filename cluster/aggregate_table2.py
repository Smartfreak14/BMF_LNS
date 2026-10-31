#!/usr/bin/env python3
"""
Agrege les CSV Table 2 (BoolMaxLS K/2, multi-seed).
Format d'entree : results_table2/<dataset>_s<N>.csv (1 CSV par seed).
Format de sortie : dataset,k,nseeds,err_min,err_median,err_max,err_std,
                   time_ms_min,time_ms_median,time_ms_max,time_to_zero_ms_median,
                   nb_zero_error,m,n,regime,tool

Usage :
    python aggregate_table2.py results_table2/ -o table2_all.csv
"""
import argparse, csv, sys, statistics, re
from pathlib import Path

PAT = re.compile(r"^(?P<dataset>.+)_s(?P<seed>\d+)\.csv$")


def collect(results_dir: Path):
    groups = {}
    for csv_path in sorted(results_dir.glob("*_s*.csv")):
        m = PAT.match(csv_path.name)
        if not m:
            continue
        dataset = m.group("dataset") + ".csv"
        with csv_path.open() as f:
            for row in csv.DictReader(f):
                if not row.get("final_errors"):
                    continue
                try:
                    err = int(row["final_errors"])
                    t   = float(row["time_ms"])
                    ttz = float(row.get("time_to_zero_ms", -1))
                except (ValueError, TypeError):
                    continue
                g = groups.setdefault(dataset, {
                    "errs": [], "times": [], "ttzs": [],
                    "k": row.get("k", ""), "m": row.get("m", ""),
                    "n": row.get("n", ""), "regime": row.get("regime", ""),
                    "tool": row.get("tool", row.get("method", "")),
                })
                g["errs"].append(err)
                g["times"].append(t)
                if ttz >= 0:
                    g["ttzs"].append(ttz)
    return groups


def main():
    p = argparse.ArgumentParser()
    p.add_argument("results_dir", type=Path)
    p.add_argument("-o", "--output", type=Path, default=Path("table2_all.csv"))
    args = p.parse_args()

    groups = collect(args.results_dir)
    if not groups:
        print(f"[error] aucun _s<N>.csv exploitable dans {args.results_dir}", file=sys.stderr)
        return 1

    header = [
        "dataset", "k", "m", "n", "regime", "tool", "nseeds",
        "err_min", "err_median", "err_max", "err_std",
        "time_ms_min", "time_ms_median", "time_ms_max", "time_ms_std",
        "time_to_zero_ms_median", "nb_zero_error",
    ]

    rows_out = []
    for d, g in sorted(groups.items()):
        errs, times, ttzs = g["errs"], g["times"], g["ttzs"]
        rows_out.append({
            "dataset": d, "k": g["k"], "m": g["m"], "n": g["n"],
            "regime": g["regime"], "tool": g["tool"],
            "nseeds": len(errs),
            "err_min": min(errs), "err_median": int(statistics.median(errs)),
            "err_max": max(errs),
            "err_std": round(statistics.pstdev(errs), 3) if len(errs) > 1 else 0,
            "time_ms_min": round(min(times), 1),
            "time_ms_median": round(statistics.median(times), 1),
            "time_ms_max": round(max(times), 1),
            "time_ms_std": round(statistics.pstdev(times), 1) if len(times) > 1 else 0,
            "time_to_zero_ms_median": round(statistics.median(ttzs), 1) if ttzs else -1,
            "nb_zero_error": sum(1 for e in errs if e == 0),
        })

    with args.output.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    print(f"[ok] {len(rows_out)} datasets ecrits dans {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
