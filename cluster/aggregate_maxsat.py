#!/usr/bin/env python3
"""
Agrege les CSV bench par-dataset (<dataset>_bench.csv) produits par
job_maxsat_array.sh en un unique CSV wide-format :

    dataset,boolmaxls_err,evalmaxsat_err,nuwls_err,satlike_err,loandra_err

Optionnellement, jointure avec un CSV boolmaxls (produit par aggregate.py)
via --with-boolmaxls <path>.

Exemple :
    python aggregate_maxsat.py results_maxsat_simpluni/ \\
        -o results_maxsat_simpluni_all.csv \\
        --with-boolmaxls results_boolmaxls_simpluni_all.csv
"""
import argparse
import csv
import sys
from pathlib import Path

SOLVERS = ["evalmaxsat", "nuwls", "satlike", "loandra"]


def collect_bench(results_dir: Path):
    rows = {}
    for csv_path in sorted(results_dir.glob("*_bench.csv")):
        dataset = csv_path.name[:-len("_bench.csv")] + ".csv"
        row = {"dataset": dataset}
        with csv_path.open() as f:
            for r in csv.DictReader(f):
                s = r["solver"]
                if s not in SOLVERS:
                    continue
                cost = r.get("best_cost") or ""
                row[f"{s}_err"] = cost if cost not in ("", "N/A") else ""
                row[f"{s}_time_ms"] = r.get("time_ms", "")
                row[f"{s}_solved"] = r.get("solved", "")
        rows[dataset] = row
    return rows


def collect_boolmaxls(path: Path):
    if not path.is_file():
        print(f"[warn] {path} absent, colonne boolmaxls_err omise", file=sys.stderr)
        return {}
    rows = {}
    with path.open() as f:
        for r in csv.DictReader(f):
            rows[r["dataset"]] = {
                "boolmaxls_err": r.get("final_errors", ""),
                "boolmaxls_time_ms": r.get("time_ms", ""),
                "boolmaxls_flips_per_sec": r.get("wls_flips_per_sec", ""),
            }
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("results_dir", type=Path)
    p.add_argument("-o", "--output", type=Path, default=Path("maxsat_baseline_all.csv"))
    p.add_argument("--with-boolmaxls", type=Path, default=None,
                   help="CSV boolmaxls agrege (aggregate.py) pour jointure.")
    args = p.parse_args()

    if not args.results_dir.is_dir():
        print(f"[error] {args.results_dir} n'est pas un dossier", file=sys.stderr)
        return 1

    bench = collect_bench(args.results_dir)
    if not bench:
        print(f"[error] aucun _bench.csv dans {args.results_dir}", file=sys.stderr)
        return 1

    opt4 = collect_boolmaxls(args.with_boolmaxls) if args.with_boolmaxls else {}

    all_datasets = sorted(set(bench) | set(opt4))
    header = ["dataset"]
    if opt4:
        header += ["boolmaxls_err", "boolmaxls_time_ms", "boolmaxls_flips_per_sec"]
    for s in SOLVERS:
        header += [f"{s}_err", f"{s}_time_ms", f"{s}_solved"]

    with args.output.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for d in all_datasets:
            row = {"dataset": d}
            row.update(bench.get(d, {}))
            row.update(opt4.get(d, {}))
            w.writerow(row)

    print(f"[ok] {len(all_datasets)} lignes ecrites dans {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
