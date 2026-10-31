#!/usr/bin/env python3
"""
Agrege les CSV Kolomvakis (dataset x r x budget x nseeds) en un CSV
au format Table 4 du papier :

    dataset,r,budget_s,err_min,err_median,err_max,err_std,nseeds

Usage :
    python aggregate_kolomvakis.py results_kolomvakis/ -o kolomvakis_all.csv
"""
import argparse, csv, sys, statistics, re
from pathlib import Path

PAT = re.compile(r"^(?P<dataset>.+)_r(?P<r>\d+)_t(?P<t>\d+)\.csv$")


def collect(results_dir: Path):
    out = []
    for csv_path in sorted(results_dir.glob("*_r*_t*.csv")):
        m = PAT.match(csv_path.name)
        if not m:
            continue
        dataset = m.group("dataset") + ".csv"
        r = int(m.group("r"))
        t = int(m.group("t"))
        errs = []
        with csv_path.open() as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    errs.append(int(row["final_errors"]))
                except (KeyError, ValueError, TypeError):
                    continue
        if not errs:
            continue
        out.append({
            "dataset": dataset,
            "r": r,
            "budget_s": t,
            "err_min": min(errs),
            "err_median": int(statistics.median(errs)),
            "err_max": max(errs),
            "err_std": round(statistics.pstdev(errs), 3) if len(errs) > 1 else 0,
            "nseeds": len(errs),
        })
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("results_dir", type=Path)
    p.add_argument("-o", "--output", type=Path, default=Path("kolomvakis_all.csv"))
    args = p.parse_args()

    rows = collect(args.results_dir)
    if not rows:
        print(f"[error] aucun CSV Kolomvakis dans {args.results_dir}", file=sys.stderr)
        return 1

    rows.sort(key=lambda x: (x["dataset"], x["r"], x["budget_s"]))
    header = list(rows[0].keys())
    with args.output.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"[ok] {len(rows)} lignes ecrites dans {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
