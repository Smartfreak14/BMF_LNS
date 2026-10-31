#!/usr/bin/env python3
"""
Depart de BoolMaxLS sur les matrices brutes, a k et floor(k/2) : erreur au budget, mediane des graines.

  GreEss          : greess_raw_k / _raw_khalf (deterministe, cellules inconnues free)
  Ours, glouton   : boolmaxls_raw_k_final / _raw_khalf_final
  Ours, GreEss    : boolmaxls_greess_raw_k / _raw_khalf, depart des facteurs de GreEss calcules dans
                    la meme tache ; l'erreur est lue au budget moins le temps de GreEss
Gras : plus petite erreur des trois colonnes, egalites toutes en gras.

Usage (depuis cluster/) :
    python make_greess_start_table.py [--root dossier_des_results] [-o table_greess_start.tex] [--budget 1800]
"""
import argparse, csv, statistics, sys
from pathlib import Path

from make_frombelow_table import at_budget, boolmaxls_runs, fmt_num, frombelow_cells, label, read_k

SETTINGS = [("raw_k", "boolmaxls_raw_k_final", "boolmaxls_greess_raw_k"),
            ("raw_khalf", "boolmaxls_raw_khalf_final", "boolmaxls_greess_raw_khalf")]


def ours(root, campaign, budget_s, greess):
    runs = {}
    for f in boolmaxls_runs(root, campaign):
        r = next(csv.DictReader(f.open()), None)
        if r is None:
            continue
        base, seed = f.stem.rsplit("_boolmaxls_s", 1)
        b = budget_s
        if greess:
            init = f.parent / f"{base}_greess_init_s{seed}" / f"{base}_greess_free.csv"
            if not init.exists():
                sys.exit(f"{init} absent : temps de GreEss inconnu")
            b -= float(next(csv.DictReader(init.open()))["total_time_s"])
        runs.setdefault(base, []).append(at_budget(f, int(r["final_errors"]), b * 1000))
    return {base: statistics.median(v) for base, v in runs.items()}


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=here)
    p.add_argument("-o", "--output", type=Path)
    p.add_argument("--budget", type=float, default=1800, help="budget par execution (s)")
    args = p.parse_args()

    ks = {"raw_k": read_k(here / "datasets_real.csv"), "raw_khalf": read_k(here / "datasets_real_khalf.csv")}
    cols = {}
    for s, greedy, greess in SETTINGS:
        for campaign in (f"greess_{s}", greedy, greess):
            if not (args.root / f"results_{campaign}").is_dir():
                sys.exit(f"results_{campaign} absent sous {args.root}")
        cols[s] = ({b: c[1][0] for b, c in frombelow_cells(args.root, f"greess_{s}", "free", False).items()
                    if c[1] is not None},
                   ours(args.root, greedy, args.budget, False),
                   ours(args.root, greess, args.budget, True))

    lines = ["\\begin{table}[h]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
             "\\begin{tabular}{l rrrr rrrr}", "\\toprule",
             " & \\multicolumn{4}{c}{raw, rank $k$} & \\multicolumn{4}{c}{raw, rank $\\lfloor k/2 \\rfloor$} \\\\",
             "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
             " & & & \\multicolumn{2}{c}{Ours, start} & & & \\multicolumn{2}{c}{Ours, start} \\\\",
             "\\textsc{Dataset} & $k$ & GreEss & greedy & GreEss "
             "& $\\lfloor k/2 \\rfloor$ & GreEss & greedy & GreEss \\\\", "\\midrule"]
    count = {s: {"<": 0, "=": 0, ">": 0} for s in cols}
    for base in sorted(ks["raw_k"], key=str.lower):
        row = [label(base)]
        for s in cols:
            vals = [c.get(base) for c in cols[s]]
            known = [v for v in vals if v is not None]
            best = min(known) if known else None
            row.append(str(ks[s][base]))
            row += ["?" if v is None else (f"\\textbf{{{fmt_num(v)}}}" if v == best else fmt_num(v)) for v in vals]
            if vals[1] is not None and vals[2] is not None:
                count[s]["<" if vals[2] < vals[1] else ("=" if vals[2] == vals[1] else ">")] += 1
        lines.append(" & ".join(row) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}",
              "\\caption{Error on the raw matrices with \\ourmethod{} started from the greedy solution "
              "or from the GreEss factors (median over the seeds, error at the 30-min limit). "
              "GreEss runs in the same 30 minutes, before \\ourmethod{}. Bold: smallest error of the row "
              "and rank, ties all bolded.}",
              "\\label{tab:greess-start}", "\\end{table}"]
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
        print(f"-> {args.output}", file=sys.stderr)
    else:
        print(text)
    for s, c in count.items():
        print(f"[{s}] depart GreEss / depart glouton : meilleur {c['<']}, egal {c['=']}, pire {c['>']}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
