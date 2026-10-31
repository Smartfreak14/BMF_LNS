#!/usr/bin/env python3
"""
Resume d'une campagne NuWLS-c (results_<RUN_NAME>/<base>_bench.csv + <base>.log).

Trois listes separees :
  - non lances      : encodage > MAX_CLAUSES (note too_large_nclauses=...)
  - sans premiere solution : lance, aucune ligne "o" (note no_first_solution;...)
  - resolus         : au moins une solution ; best_cost = erreur finale.

Depart glouton (log "[init] greedy errors=E") : verifie first_cost == E et best <= E.

Usage :
    python analyze_nuwls_runs.py results_nuwls_greedy_raw_k/ datasets_real.csv
"""
import argparse, csv, re, sys
from pathlib import Path


def read_k(path):
    rows = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")][1:]
    return {r.split(",")[0].strip(): int(r.split(",")[1]) for r in rows}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("results_dir", type=Path)
    p.add_argument("datasets", type=Path)
    args = p.parse_args()

    k_list = read_k(args.datasets)
    missing, not_launched, no_first, solved = [], [], [], []
    problems = 0

    print(f"{'dataset':24s}{'k':>6s}{'glouton':>10s}{'1re sol':>10s}{'t 1re (s)':>11s}{'best':>10s}  statut")
    for ds in sorted(k_list, key=str.lower):
        base = ds[:-4]
        bench = args.results_dir / f"{base}_bench.csv"
        if not bench.exists():
            missing.append(ds)
            print(f"{ds:24s}{k_list[ds]:>6d}{'':>41s}  MANQUANT")
            continue
        with bench.open() as f:
            r = next((x for x in csv.DictReader(f) if x["solver"] == "nuwls"), None)
        if r is None:
            missing.append(ds)
            continue
        log = args.results_dir / f"{base}.log"
        m = re.search(r"\[init\] greedy errors=(\d+)", log.read_text(errors="ignore")) if log.exists() else None
        greedy = m.group(1) if m else ""
        note = r.get("note", "")
        first = r.get("first_cost") or ""
        msg = ""
        if note.startswith("too_large"):
            not_launched.append(f"{base} ({note.split('=')[1]} clauses)")
            status = "non lance (> seuil)"
        elif note.startswith("no_first_solution"):
            no_first.append(f"{base} ({note.split(';', 1)[1]})")
            status = "sans premiere solution"
        else:
            solved.append(ds)
            status = "ok"
            if greedy and first != greedy:
                msg += f" ; 1re sol {first} != glouton {greedy}"
            if greedy and int(r["best_cost"]) > int(greedy):
                msg += " ; best > glouton"
        problems += bool(msg)
        print(f"{ds:24s}{k_list[ds]:>6d}{greedy:>10s}{first:>10s}{r.get('first_time_s') or '':>11s}"
              f"{r['best_cost']:>10s}  {status}{msg}")

    print(f"\nresolus : {len(solved)}/{len(k_list)}")
    print(f"non lances (encodage > seuil) : {len(not_launched)}")
    for x in not_launched:
        print(f"  {x}")
    print(f"lances sans premiere solution : {len(no_first)}")
    for x in no_first:
        print(f"  {x}")
    if missing:
        print(f"manquants : {len(missing)} -> {', '.join(missing)}")
    problems += len(missing)

    print("\n[ok] campagne complete et coherente" if problems == 0 else f"\n[warn] {problems} probleme(s)")
    return 0 if problems == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
