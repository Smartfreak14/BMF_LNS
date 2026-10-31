#!/usr/bin/env python3
"""
Tableau LaTeX des jeux de donnees bruts : m, n, densite, 1 connus, cellules inconnues, k, min(m,n).

Lecture identique a CSVMatrixLoader : un champ vide ou "?" est une cellule inconnue (y compris
le dernier champ d'une ligne terminee par une virgule) ; densite = 1 connus / (m n).
k = colonne k de datasets_real.csv (rang trouve par OptiBlock sur la matrice originale,
Avellaneda et Villemaire, AAAI 2024, tableau 2).

Usage (depuis cluster/) :
    python make_dataset_table.py [-o table_datasets.tex] [--data ../data/Real]
"""
import argparse, sys
from pathlib import Path

LABELS = {"americas_small": "Americas\\_small", "chess_krvskp": "Chess",
          "studentPerformance": "StudentPerf.", "ThoraricSurgery": "Thoracic"}


def stats(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.strip():
            rows.append(line.rstrip("\r").split(","))
    m, n = len(rows), max(len(r) for r in rows)
    ones = unknown = 0
    for r in rows:
        for t in r:
            t = t.strip()
            if t in ("", "?"):
                unknown += 1
            elif float(t) != 0:
                ones += 1
    return m, n, ones, unknown


def num(x):
    return f"{x:,}".replace(",", "\\,")


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=here.parent / "data" / "Real")
    p.add_argument("--datasets", type=Path, default=here / "datasets_real.csv")
    p.add_argument("-o", "--output", type=Path)
    args = p.parse_args()

    rows = [l for l in args.datasets.read_text().splitlines() if l.strip() and not l.startswith("#")][1:]
    lines = ["\\begin{table}[t]", "\\centering", "\\scriptsize",
             "\\setlength{\\tabcolsep}{3pt}",
             "\\begin{tabular}{l rrrrrrr}", "\\toprule",
             "\\textsc{Dataset} & $m$ & $n$ & density & known 1s & unknown & $k$ & $\\min(m,n)$ \\\\",
             "\\midrule"]
    for r in sorted(rows, key=str.lower):
        f, k = r.split(",")[0].strip(), int(r.split(",")[1])
        base = f[:-4]
        m, n, ones, unknown = stats(args.data / f)
        name = LABELS.get(base, base[0].upper() + base[1:])
        lines.append(f"{name} & {num(m)} & {num(n)} & {ones / (m * n):.3f} & {num(ones)} & "
                     f"{num(unknown) if unknown else '--'} & {num(k)} & {num(min(m, n))} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}",
              "\\caption{The 31 raw benchmark matrices: $m$ rows, $n$ columns, density = known 1s$/(mn)$, "
              "number of known 1s, number of unknown cells (``--'' = none), rank $k$ used in the raw "
              "experiments (rank of the factorization found by OptiBlock on the original matrix within "
              "3\\,h, Avellaneda and Villemaire 2024, Table~2) and $\\min(m,n)$.}",
              "\\label{tab:datasets}", "\\end{table}"]
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
        print(f"-> {args.output}", file=sys.stderr)
    else:
        print(text)


if __name__ == "__main__":
    main()
