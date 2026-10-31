#!/usr/bin/env python3
"""
Tableau LaTeX de l'ablation sur les matrices brutes, a k et floor(k/2), au format de tab:ablation :
erreur au budget (d'apres la trace) des cinq variantes emboitees GI (greedy), no-w (no-weighting),
+wt (wls-no-double-no-lns), +DF (no-lns) et full, mediane des graines (ablation_raw_k / _raw_khalf).
Gras : plus petite erreur de la ligne pour chaque rang, egalites toutes en gras ; derniere ligne :
nombre de matrices ou chaque variante est en gras.
Sur stderr : executions manquantes, tuees (verif = -1), verifiees en desaccord, ou finies apres le budget.

Usage (depuis cluster/) :
    python make_ablation_table.py [--root dossier_des_results] [-o table_ablation.tex] [--budget 1800]
"""
import argparse, csv, re, statistics, sys
from pathlib import Path

from make_frombelow_table import at_budget, fmt_num, label, read_k

VARIANTS = [("greedy", "GI"), ("no-weighting", "no-w"), ("wls-no-double-no-lns", "+wt"),
            ("no-lns", "+DF"), ("full", "full")]
PAT = re.compile(r"^(?P<base>.+)_(?P<mode>" + "|".join(re.escape(v) for v, _ in VARIANTS) + r")_s(?P<seed>\d+)$")
SETTINGS = [("ablation_raw_k", "datasets_real.csv"), ("ablation_raw_khalf", "datasets_real_khalf.csv")]


def runs(d, budget_s, issues):
    out = {}
    for f in sorted(d.glob("*.csv")):
        m = PAT.match(f.stem)
        r = next(csv.DictReader(f.open()), None) if m else None
        if r is None:
            continue
        name = f"{d.name[len('results_'):]} {f.stem}"
        if r["verif_errors"] == "-1":
            issues.append(f"{name} : non verifiee (tuee ou verification interrompue), meilleur {r['final_errors']} "
                          f"ecrit a {float(r['time_ms']) / 1e3:.0f}s")
        elif r["verif_errors"] != r["final_errors"]:
            issues.append(f"{name} : verif={r['verif_errors']} final={r['final_errors']}")
        elif float(r["time_ms"]) > (budget_s + 5) * 1000:
            issues.append(f"{name} : fin a {float(r['time_ms']) / 1e3:.0f}s")
        out.setdefault((m["base"], m["mode"]), []).append(at_budget(f, int(r["final_errors"]), budget_s * 1000))
    return {key: statistics.median(v) for key, v in out.items()}, max((len(v) for v in out.values()), default=0)


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=here)
    p.add_argument("-o", "--output", type=Path)
    p.add_argument("--budget", type=float, default=1800, help="budget par execution (s)")
    args = p.parse_args()

    issues, data, ks, nseeds = [], [], [], 0
    for campaign, lst in SETTINGS:
        d = args.root / f"results_{campaign}"
        if not d.is_dir():
            sys.exit(f"results_{campaign} absent sous {args.root}")
        vals, n = runs(d, args.budget, issues)
        data.append(vals)
        ks.append(read_k(here / lst))
        nseeds = max(nseeds, n)

    heads = " & ".join(n for _, n in VARIANTS)
    lines = ["\\begin{table*}[t]", "\\centering",
             "\\setlength{\\tabcolsep}{3.2pt}\\renewcommand{\\arraystretch}{0.92}", "\\scriptsize",
             "\\begin{tabular}{l r rrrrr r rrrrr}", "\\toprule",
             " & \\multicolumn{6}{c}{raw, rank $k$} & \\multicolumn{6}{c}{raw, rank $\\lfloor k/2 \\rfloor$} \\\\",
             "\\cmidrule(lr){2-7}\\cmidrule(lr){8-13}",
             f"\\textsc{{Dataset}} & $k$ & {heads} & $\\lfloor k/2 \\rfloor$ & {heads} \\\\", "\\midrule"]
    best = [[0] * len(VARIANTS) for _ in SETTINGS]
    for base in sorted(ks[0], key=str.lower):
        row = [label(base)]
        for s, vals in enumerate(data):
            v = [vals.get((base, mode)) for mode, _ in VARIANTS]
            for mode, _ in VARIANTS:
                if (base, mode) not in vals:
                    issues.append(f"{SETTINGS[s][0]} {base}_{mode} : manquante")
            known = [x for x in v if x is not None]
            low = min(known) if known else None
            row.append(str(ks[s][base]))
            for i, x in enumerate(v):
                if x is not None and x == low:
                    best[s][i] += 1
                row.append("?" if x is None else (f"\\textbf{{{fmt_num(x)}}}" if x == low else fmt_num(x)))
        lines.append(" & ".join(row) + " \\\\")
    lines += ["\\midrule",
              "\\emph{best on} & & " + " & ".join(map(str, best[0])) + " & & " + " & ".join(map(str, best[1])) + " \\\\",
              "\\bottomrule", "\\end{tabular}",
              "\\caption{Component ablation on the raw matrices at rank $k$ and $\\lfloor k/2 \\rfloor$: error at the "
              "30-minute limit of five nested variants (GI $\\subset$ no-w $\\subset$ +wt $\\subset$ +DF $\\subset$ "
              "full), " + (f"median over {nseeds} seeds" if nseeds > 1 else "one run each") + ". "
              "Bold: smallest error per row and rank, ties all bolded. "
              "Last row: number of datasets on which each variant is in bold.}",
              "\\label{tab:ablation-raw}", "\\end{table*}"]
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
        print(f"-> {args.output}", file=sys.stderr)
    else:
        print(text)
    for i in issues:
        print(f"[warn] {i}", file=sys.stderr)


if __name__ == "__main__":
    main()
