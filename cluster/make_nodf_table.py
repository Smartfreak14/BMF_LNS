#!/usr/bin/env python3
"""
Methode complete contre methode complete sans double flip, sur les 31 matrices brutes a k et
floor(k/2) : erreur au budget (d'apres la trace), mediane des graines.

  full : boolmaxls_raw_k_final / _raw_khalf_final (executions de <campagne>_redo comprises)
  -DF  : boolmaxls_nodf_raw_k / _raw_khalf (--ablation no-double-flip, memes reglages ; l'appel LNS
         en cours est coupe a la fin du budget, ce qui ne change pas l'erreur a 30 min)
Gras : plus petite des deux medianes, egalites toutes en gras. Sur stderr, par rang : nombre de
matrices ou retirer le double flip degrade, ameliore ou ne change pas la mediane, dont les ecarts
ou toutes les graines d'une variante font mieux que toutes celles de l'autre ; puis les executions
manquantes, non verifiees, verifiees en desaccord ou finies apres le budget.

Usage (depuis cluster/) :
    python make_nodf_table.py [--root dossier_des_results] [-o table_nodf.tex] [--budget 1800]
"""
import argparse, csv, statistics, sys
from pathlib import Path

from make_frombelow_table import at_budget, boolmaxls_runs, fmt_num, label, read_k

SETTINGS = [("raw_k", "boolmaxls_raw_k_final", "boolmaxls_nodf_raw_k", "datasets_real.csv"),
            ("raw_khalf", "boolmaxls_raw_khalf_final", "boolmaxls_nodf_raw_khalf", "datasets_real_khalf.csv")]


def errors(root, campaign, budget_s, issues):
    """{matrice: [erreur au budget par graine]}."""
    runs = {}
    for f in boolmaxls_runs(root, campaign):
        r = next(csv.DictReader(f.open()), None)
        if r is None:
            continue
        name = f"{f.parent.name[len('results_'):]} {f.stem}"
        if r["verif_errors"] == "-1":
            issues.append(f"{name} : non verifiee (tuee ou verification interrompue), meilleur {r['final_errors']} "
                          f"ecrit a {float(r['time_ms']) / 1e3:.0f}s")
        elif r["verif_errors"] != r["final_errors"]:
            issues.append(f"{name} : verif={r['verif_errors']} final={r['final_errors']}")
        elif float(r["time_ms"]) > (budget_s + 5) * 1000:
            issues.append(f"{name} : fin a {float(r['time_ms']) / 1e3:.0f}s")
        base = f.stem.rsplit("_boolmaxls_s", 1)[0]
        runs.setdefault(base, []).append(at_budget(f, int(r["final_errors"]), budget_s * 1000))
    return runs


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=here)
    p.add_argument("-o", "--output", type=Path)
    p.add_argument("--budget", type=float, default=1800, help="budget par execution (s)")
    args = p.parse_args()

    issues, data, ks, nseeds = [], [], [], 0
    for s, full, nodf, lst in SETTINGS:
        for campaign in (full, nodf):
            if not (args.root / f"results_{campaign}").is_dir():
                sys.exit(f"results_{campaign} absent sous {args.root}")
        runs = (errors(args.root, full, args.budget, issues), errors(args.root, nodf, args.budget, issues))
        data.append(runs)
        ks.append(read_k(here / lst))
        nseeds = max([nseeds] + [len(v) for r in runs for v in r.values()])

    lines = ["\\begin{table}[h]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
             "\\begin{tabular}{l rrr rrr}", "\\toprule",
             " & \\multicolumn{3}{c}{raw, rank $k$} & \\multicolumn{3}{c}{raw, rank $\\lfloor k/2 \\rfloor$} \\\\",
             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
             "\\textsc{Dataset} & $k$ & full & $-$DF & $\\lfloor k/2 \\rfloor$ & full & $-$DF \\\\", "\\midrule"]
    count = [{"worse": [], "better": [], "equal": [], "clear_worse": [], "clear_better": []} for _ in SETTINGS]
    best = [[0, 0] for _ in SETTINGS]
    for base in sorted(ks[0], key=str.lower):
        row = [label(base)]
        for s, (full, nodf) in enumerate(data):
            a, b = full.get(base), nodf.get(base)
            for v, campaign in ((a, SETTINGS[s][1]), (b, SETTINGS[s][2])):
                if v is None:
                    issues.append(f"{campaign} {base} : manquante")
                elif len(v) < nseeds:
                    issues.append(f"{campaign} {base} : {len(v)}/{nseeds} graines")
            ma = statistics.median(a) if a else None
            mb = statistics.median(b) if b else None
            low = min(x for x in (ma, mb) if x is not None) if (a or b) else None
            row.append(str(ks[s][base]))
            for i, x in enumerate((ma, mb)):
                if x is not None and x == low:
                    best[s][i] += 1
                row.append("?" if x is None else (f"\\textbf{{{fmt_num(x)}}}" if x == low else fmt_num(x)))
            if a and b:
                c = count[s]
                c["worse" if mb > ma else "better" if mb < ma else "equal"].append(label(base))
                if min(b) > max(a):
                    c["clear_worse"].append(label(base))
                elif max(b) < min(a):
                    c["clear_better"].append(label(base))
        lines.append(" & ".join(row) + " \\\\")
    lines += ["\\midrule",
              f"\\emph{{best on}} & & {best[0][0]} & {best[0][1]} & & {best[1][0]} & {best[1][1]} \\\\",
              "\\bottomrule", "\\end{tabular}",
              "\\caption{Error at the 30-minute limit of the full method and of the full method without the "
              "double flip ($-$DF) on the raw matrices, "
              + (f"median over {nseeds} seeds" if nseeds > 1 else "one run each") + ". "
              "Bold: smaller error, ties both bolded.}",
              "\\label{tab:no-double-flip}", "\\end{table}"]
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
        print(f"-> {args.output}", file=sys.stderr)
    else:
        print(text)

    total = {key: sum(len(c[key]) for c in count) for key in count[0]}
    for (s, *_), c in zip(SETTINGS, count):
        n = sum(len(c[key]) for key in ("worse", "better", "equal"))
        print(f"[{s}] sans double flip, sur {n} matrices : pire {len(c['worse'])} "
              f"(toutes graines : {len(c['clear_worse'])}), meilleur {len(c['better'])} "
              f"(toutes graines : {len(c['clear_better'])}), egal {len(c['equal'])}", file=sys.stderr)
        for key in ("worse", "better"):
            if c[key]:
                print(f"    {'pire' if key == 'worse' else 'meilleur'} : {', '.join(c[key])}", file=sys.stderr)
    n = total["worse"] + total["better"] + total["equal"]
    print(f"[phrase] Removing the double flip worsens the median error on {total['worse']} and improves it "
          f"on {total['better']} of the {n} instances ({total['equal']} ties).", file=sys.stderr)
    for i in issues:
        print(f"[warn] {i}", file=sys.stderr)


if __name__ == "__main__":
    main()
