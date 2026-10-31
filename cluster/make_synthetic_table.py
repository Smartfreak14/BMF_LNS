#!/usr/bin/env python3
"""
Tableau LaTeX des 12 instances synthetiques (datasets_synthetic.csv) : temps pour atteindre une
factorisation exacte au rang k plante.

  Ours                 : boolmaxls_synthetic, mediane des graines [min-max] (time_to_zero_ms)
  NuWLS                : nuwls_greedy_synthetic, NuWLS-c demarre du glouton (premiere erreur 0 de la trace)
  GreConD, GreEss, MDL : <algo>_synthetic, deterministes (une execution)
Une graine sans factorisation exacte dans le budget compte comme un temps infini dans la mediane ;
le nombre de graines exactes est alors indique. "--" = mediane infinie, "?" = campagne absente.
Gras : plus petit temps median (tel qu'affiche), egalites toutes en gras.

Usage (depuis cluster/) :
    python make_synthetic_table.py [--root dossier_des_results] [-o table_synthetic.tex] [--budget 1800]
"""
import argparse, csv, math, re, sys
from pathlib import Path

from make_frombelow_table import fmt_time, shown, read_k

FROMBELOW = [("grecond", "GreConD"), ("greess", "GreEss"), ("mdlgrecond", "MDL")]


def boolmaxls_times(d, budget_s):
    runs = {}
    for f in d.glob("*_boolmaxls_s*.csv"):
        r = next(csv.DictReader(f.open()), None)
        if r:
            t = float(r["time_to_zero_ms"]) / 1e3
            ok = int(r["final_errors"]) == 0 and 0 <= t <= budget_s
            runs.setdefault(f.name.rsplit("_boolmaxls_s", 1)[0], []).append(t if ok else math.inf)
    return runs


def nuwls_times(d, budget_s):
    runs = {}
    for f in d.glob("*_bench.csv"):
        base = re.sub(r"_s\d+$", "", f.name[:-len("_bench.csv")])
        tr = f.with_name(f.name[:-len("_bench.csv")] + "_nuwls.trace")
        pts = [l.split() for l in tr.read_text().splitlines()[1:]] if tr.exists() else []
        t = next((float(p[0]) / 1e3 for p in pts if len(p) == 2 and int(p[1]) == 0), math.inf)
        runs.setdefault(base, []).append(t if t <= budget_s else math.inf)
    return runs


def frombelow_times(d, algo):
    runs = {}
    for f in d.glob(f"*_{algo}_free.csv"):
        r = next(csv.DictReader(f.open()), None)
        if r:
            ok = r["error"] == "0" and not r["stop"].startswith("timeout")
            runs[r["dataset"]] = [float(r["time_s"]) if ok else math.inf]
    return runs


def cell(ts, spread):
    """(texte, cle de comparaison) ; mediane basse avec un nombre pair de graines."""
    if ts is None:
        return ("?", None)
    ts = sorted(ts)
    m = ts[(len(ts) - 1) // 2]
    if math.isinf(m):
        return ("--", None)
    text = fmt_time(m)
    ok = [t for t in ts if not math.isinf(t)]
    if spread and len(ts) > 1:
        text += f" [{fmt_time(ok[0])}--{fmt_time(ok[-1])}]"
    if len(ok) < len(ts):
        text += f" ({len(ok)}/{len(ts)})"
    return (text, shown(m))


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=here)
    p.add_argument("-o", "--output", type=Path)
    p.add_argument("--budget", type=float, default=1800, help="budget par execution (s)")
    args = p.parse_args()

    ks = read_k(here / "datasets_synthetic.csv")
    cols = [("Ours", "boolmaxls_synthetic", lambda d: boolmaxls_times(d, args.budget), True),
            ("NuWLS", "nuwls_greedy_synthetic", lambda d: nuwls_times(d, args.budget), False)]
    cols += [(name, f"{a}_synthetic", lambda d, a=a: frombelow_times(d, a), False) for a, name in FROMBELOW]
    data, nseeds = {}, {}
    for name, campaign, read, _ in cols:
        d = args.root / f"results_{campaign}"
        data[name] = read(d) if d.is_dir() else None
        if data[name] is not None:
            nseeds[name] = max((len(v) for v in data[name].values()), default=0)
        print(f"[source] {name:7s} : " + (f"results_{campaign} ({len(data[name])}/{len(ks)} instances)"
                                          if data[name] is not None else f"results_{campaign} absent (\"?\")"),
              file=sys.stderr)

    lines = ["\\begin{table}[h]", "\\centering", "\\small",
             f"\\begin{{tabular}}{{rr {'r' * len(cols)}}}", "\\toprule",
             "density & $k$ & " + " & ".join(n for n, *_ in cols) + " \\\\", "\\midrule"]
    for base in sorted(ks, key=lambda b: (ks[b], b)):
        m = re.match(r"^d_([\d.]+)_k_(\d+)$", base)
        cells = [cell(None if data[n] is None else data[n].get(base), spread) for n, _, _, spread in cols]
        keys = [c[1] for c in cells if c[1] is not None]
        best = min(keys) if keys else None
        txt = [f"\\textbf{{{t}}}" if k is not None and k == best else t for t, k in cells]
        lines.append(f"{m.group(1) if m else base} & {ks[base]} & " + " & ".join(txt) + " \\\\")
    seeds = nseeds.get("Ours", 0)
    lines += ["\\bottomrule", "\\end{tabular}",
              f"\\caption{{Time to an exact factorization at the planted rank $k$ on the synthetic "
              f"$1000 \\times 500$ matrices (30-min limit). \\emph{{Ours}}: median over {seeds} seeds "
              "[min--max]; \\emph{NuWLS} = NuWLS-c started from the greedy solution, median over the same "
              "seeds, with the number of seeds reaching an exact factorization when not all do; GreConD, "
              "GreEss and MDL (MDLGreConD) are deterministic (one run). ``--'' = no exact factorization "
              "within the limit for at least half of the runs. Bold: shortest median time.}",
              "\\label{tab:synthetic}", "\\end{table}"]
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
        print(f"-> {args.output}", file=sys.stderr)
    else:
        print(text)


if __name__ == "__main__":
    main()
