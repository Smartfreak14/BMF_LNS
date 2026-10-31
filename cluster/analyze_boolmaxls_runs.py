#!/usr/bin/env python3
"""
Controle et resume d'une campagne BoolMaxLS (results_<RUN_NAME>/<base>_boolmaxls_s<N>.csv).

  1. Completude : un CSV par (dataset, graine).
  2. k du CSV == k de la liste de datasets.
  3. Verification : verif_errors == final_errors (verif = -1 -> verificateur non execute).
  4. From-below : aucun "from-below VIOLATED" dans les .log.
  5. Budget : time_ms <= budget + 5 s.

Usage :
    python analyze_boolmaxls_runs.py results_boolmaxls_raw_k/ datasets_real.csv [--nseeds 1] [--budget 1800]
"""
import argparse, csv, statistics, sys
from pathlib import Path


def read_k(path):
    rows = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")][1:]
    return {r.split(",")[0].strip(): int(r.split(",")[1]) for r in rows}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("results_dir", type=Path)
    p.add_argument("datasets", type=Path)
    p.add_argument("--nseeds", type=int, default=1)
    p.add_argument("--budget", type=float, default=1800)
    args = p.parse_args()

    k_list = read_k(args.datasets)
    problems = 0
    print(f"{'dataset':24s}{'k':>6s}{'init':>9s}{'final':>9s}{'time_s':>9s}{'verif':>8s}")
    for ds in sorted(k_list, key=str.lower):
        base = ds[:-4]
        finals = []
        for s in range(1, args.nseeds + 1):
            f = args.results_dir / f"{base}_boolmaxls_s{s}.csv"
            if not f.exists():
                print(f"{ds:24s}  manquant s{s}")
                problems += 1
                continue
            with f.open() as fh:
                r = next(csv.DictReader(fh))
            msgs = []
            if int(r["k"]) != k_list[ds]:
                msgs.append(f"k={r['k']} attendu {k_list[ds]}")
            if r["verif_errors"] != r["final_errors"]:
                msgs.append(f"verif={r['verif_errors']} final={r['final_errors']}")
            if float(r["time_ms"]) > (args.budget + 5) * 1000:
                msgs.append("depassement budget")
            log = f.with_suffix(".log")
            if log.exists() and "from-below VIOLATED" in log.read_text(errors="ignore"):
                msgs.append("FROM-BELOW VIOLE")
            problems += len(msgs)
            finals.append(int(r["final_errors"]))
            print(f"{ds:24s}{int(r['k']):>6d}{r['init_errors']:>9s}{r['final_errors']:>9s}"
                  f"{float(r['time_ms']) / 1000:>9.1f}{r['verif_errors']:>8s}  s{s} {' ; '.join(msgs)}")
        if len(finals) > 1:
            print(f"{'':24s}  mediane={statistics.median(finals):g} min={min(finals)}")

    print("\n[ok] aucun probleme detecte" if problems == 0 else f"\n[warn] {problems} probleme(s)")
    return 0 if problems == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
