#!/usr/bin/env python3
"""
Controle et resume de la campagne d'ablation k/2 (results_ablation_real_khalf/).

  1. Completude : 31 datasets x 5 modes x NSEEDS CSV attendus.
  2. k du CSV de sortie == k de datasets_real_khalf.csv == floor(k_real/2).
  3. From-below : aucun "from-below VIOLATED" dans les .log, verif_errors == final_errors.
  4. Tableau : mediane des erreurs par mode et nombre de datasets ou chaque mode est le meilleur.

Usage :
    python analyze_ablation_khalf.py results_ablation_real_khalf/ [--nseeds 5]
"""
import argparse, csv, re, statistics, sys
from pathlib import Path

MODES = ["greedy", "no-weighting", "wls-no-double-no-lns", "no-lns", "full"]
PAT = re.compile(r"^(?P<ds>.+?)_(?P<mode>wls-no-double-no-lns|no-weighting|greedy|no-lns|full)_s(?P<seed>\d+)$")


def read_k(path):
    rows = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")][1:]
    return {r.split(",")[0].strip(): int(r.split(",")[1]) for r in rows}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("results_dir", type=Path)
    p.add_argument("--nseeds", type=int, default=1)
    args = p.parse_args()

    here = Path(__file__).resolve().parent
    k_real = read_k(here / "datasets_real.csv")
    k_half = read_k(here / "datasets_real_khalf.csv")

    runs = {}
    for f in sorted(args.results_dir.glob("*.csv")):
        m = PAT.match(f.stem)
        if not m:
            continue
        with f.open() as fh:
            r = next(csv.DictReader(fh), None)
        if r:
            runs[(m["ds"] + ".csv", m["mode"], int(m["seed"]))] = r

    problems = 0

    missing = [(ds, mo, s) for ds in k_half for mo in MODES for s in range(1, args.nseeds + 1)
               if (ds, mo, s) not in runs]
    expected = len(k_half) * len(MODES) * args.nseeds
    print(f"[completude] {expected - len(missing)}/{expected} runs")
    for ds, mo, s in missing[:20]:
        print(f"  manquant : {ds} {mo} s{s}")
    problems += len(missing)

    for (ds, mo, s), r in runs.items():
        k = int(r["k"])
        if k != k_half.get(ds) or k != k_real[ds] // 2:
            print(f"[k] {ds} {mo} s{s} : k={k} attendu {k_real[ds] // 2}")
            problems += 1
        if r["verif_errors"] != r["final_errors"]:
            print(f"[verif] {ds} {mo} s{s} : verif={r['verif_errors']} final={r['final_errors']}")
            problems += 1

    violated = [f.name for f in args.results_dir.glob("*.log") if "from-below VIOLATED" in f.read_text(errors="ignore")]
    print(f"[from-below] {len(violated)} run(s) avec 0-entries couvertes")
    for name in violated:
        print(f"  {name}")
    problems += len(violated)

    print()
    print(f"{'dataset':24s}{'k':>6s}" + "".join(f"{mo[:12]:>14s}" for mo in MODES))
    wins = {mo: 0 for mo in MODES}
    for ds in sorted(k_half, key=str.lower):
        med = {}
        for mo in MODES:
            errs = [int(runs[(ds, mo, s)]["final_errors"]) for s in range(1, args.nseeds + 1) if (ds, mo, s) in runs]
            if errs:
                med[mo] = statistics.median(errs)
        if med:
            best = min(med.values())
            for mo, v in med.items():
                wins[mo] += v == best
        print(f"{ds:24s}{k_half[ds]:>6d}" + "".join(f"{med.get(mo, float('nan')):>14g}" for mo in MODES))
    print()
    print("meilleur (ex aequo inclus) : " + ", ".join(f"{mo}={wins[mo]}" for mo in MODES))

    print("\n[ok] aucun probleme detecte" if problems == 0 else f"\n[warn] {problems} probleme(s)")
    return 0 if problems == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
