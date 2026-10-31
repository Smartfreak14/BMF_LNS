#!/usr/bin/env python3
"""
Resume d'une campagne GreConD / GreEss / MDLGreConD (job_frombelow_dataset.sh).

Pour chaque matrice et chaque traitement des cellules inconnues (free, zero, one) :
facteurs, erreur, temps, et controles
  - auto-controle du programme (produit A o B recalcule) = ok
  - verif.py : meme nombre de facteurs, meme erreur, 0 zero couvert (from-below)
Listes separees : arrets par le budget, facteurs vides (MDLGreConD), manquants.
Pour les campagnes jusqu'a l'exact (simpluni, simplexi), "exact" = erreur nulle.

Usage :
    python analyze_frombelow_runs.py results_grecond_raw_k/ datasets_real.csv [-o grecond_raw_k_all.csv]
"""
import argparse, csv, sys
from pathlib import Path

MODES = ("free", "zero", "one")


def read_k(path):
    rows = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")][1:]
    return {r.split(",")[0].strip()[:-4]: int(r.split(",")[1]) for r in rows}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("results_dir", type=Path)
    p.add_argument("datasets", type=Path)
    p.add_argument("-o", "--output", type=Path)
    args = p.parse_args()

    k_list = read_k(args.datasets)
    missing, problems, timeouts, empty, out_rows = [], [], [], [], []

    print(f"{'dataset':22s}{'inconnues':>10s}{'k liste':>9s}{'facteurs':>10s}{'erreur':>9s}{'t (s)':>10s}  arret")
    for base in sorted(k_list, key=str.lower):
        found = False
        for mode in MODES:
            f = next(args.results_dir.glob(f"{base}_*_{mode}.csv"), None)
            if f is None or f.stem.count("_") != base.count("_") + 2:
                continue
            found = True
            r = next(csv.DictReader(f.open()), None)
            if r is None:
                problems.append(f"{base} {mode} : resultat vide")
                continue
            vpath = f.with_name(f.stem + "_verif.csv")
            v = next(csv.DictReader(vpath.open()), None) if vpath.exists() else None
            msg = []
            if r["self_check"] != "ok":
                msg.append("auto-controle")
            if int(r["covered_zeros"]) != 0:
                msg.append("zeros couverts")
            if int(r["factors"]) > 0:
                if v is None:
                    msg.append("verif absente")
                elif (v["verif_k"] != r["factors"] or v["error"] != r["error"]
                      or v["covered_zeros"] != "0"):
                    msg.append(f"verif {v['verif_k']}/{v['error']}/{v['covered_zeros']}")
            if msg:
                problems.append(f"{base} {mode} : {', '.join(msg)}")
            if r["stop"].startswith("timeout"):
                timeouts.append(f"{base} {mode} : {r['stop']} apres {r['factors']} facteurs")
            if r["nonempty_factors"] != r["factors"]:
                empty.append(f"{base} {mode} : {r['nonempty_factors']} non vide(s) sur {r['factors']}")
            out_rows.append({**r, "verif_error": v["error"] if v else "",
                             "verif_covered_zeros": v["covered_zeros"] if v else ""})
            stop = r["stop"] + (" (exact)" if r["error"] == "0" else "")
            print(f"{base:22s}{mode:>10s}{k_list[base]:>9d}{r['factors']:>10s}{r['error']:>9s}"
                  f"{r['time_s']:>10s}  {stop}")
        if not found:
            missing.append(base)
            print(f"{base:22s}  MANQUANT")

    print(f"\narrets par le budget : {len(timeouts)}")
    for t in timeouts:
        print(f"  {t}")
    # MDLGreConD peut retenir le concept vide (extension la moins couteuse en bits).
    print(f"facteurs vides (extent ou intent vide) : {len(empty)}")
    for t in empty:
        print(f"  {t}")
    if missing:
        print(f"manquants : {len(missing)} -> {', '.join(missing)}")
    for pb in problems:
        print(f"  [pb] {pb}")
    if args.output and out_rows:
        with args.output.open("w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(out_rows[0]))
            w.writeheader()
            w.writerows(out_rows)
        print(f"-> {args.output}")
    n_bad = len(missing) + len(problems)
    print("\n[ok] campagne complete, toutes les sorties verifiees" if n_bad == 0 else f"\n[warn] {n_bad} probleme(s)")
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
