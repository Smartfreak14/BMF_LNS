#!/usr/bin/env python3
"""
Courbes erreur-temps a partir des fichiers .trace des campagnes de launch_rerun.sh.

  raw_k / raw_khalf : Ours (boolmaxls_raw_k[half]), NuWLS (nuwls_greedy_raw_k[half]),
                      GreConD, GreEss, MDL (<algo>_raw_k[half], cellules inconnues free)
  ablation          : les 5 variantes de ablation_raw_khalf (GI, no-w, +wt, +DF, full)
  ablation_k        : les memes sur ablation_raw_k

Une page par matrice (erreur en fonction du temps, temps en echelle log, palier apres chaque
amelioration, prolonge jusqu'a la fin de l'execution) : figures/<reglage>.pdf (multipage) et
figures/<reglage>/<matrice>.pdf. Toutes les graines presentes sont tracees.

Usage (depuis cluster/, module load scipy-stack sur Narval) :
    python plot_traces.py raw_k [--root .] [--out figures] [--datasets zoo,votes]
"""
import argparse, csv, sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3de", "#fcfcfb"
# Couleurs categorielles dans l'ordre fixe de la palette validee, plus un style de trait.
METHODS = [("Ours", "#2a78d6", "-", "o"), ("NuWLS", "#eb6834", "--", "s"),
           ("GreConD", "#1baf7a", ":", "^"), ("GreEss", "#eda100", "-.", "D"),
           ("MDL", "#e87ba4", (0, (5, 1, 1, 1)), "v")]
ABLATION = [("full", "full", "#2a78d6", "-", "o"), ("greedy", "GI", "#eb6834", "--", "s"),
            ("no-weighting", "no-w", "#1baf7a", ":", "^"),
            ("wls-no-double-no-lns", "+wt", "#eda100", "-.", "D"),
            ("no-lns", "+DF", "#e87ba4", (0, (5, 1, 1, 1)), "v")]
LABELS = {"americas_small": "Americas_small", "chess_krvskp": "Chess",
          "studentPerformance": "StudentPerf.", "ThoraricSurgery": "Thoracic"}
T_MIN = 1e-4  # s, plancher de l'axe log


def read_list(path):
    rows = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")][1:]
    return [(r.split(",")[0].strip()[:-4], int(r.split(",")[1])) for r in rows]


def read_trace(path):
    """[(temps_s, erreur)] ; formats 'time_ms err' (BoolMaxLS, NuWLS) et CSV des baselines."""
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    if not lines:
        return []
    if lines[0].startswith("factors,"):
        return [(float(r["time_s"]), int(r["uncovered_ones"])) for r in csv.DictReader(lines)]
    return [(float(a) / 1e3, int(b)) for a, b in (l.split()[:2] for l in lines[1:])]


def end_time(csv_path, column, scale):
    if not csv_path.exists():
        return None
    r = next(csv.DictReader(csv_path.open()), None)
    try:
        return float(r[column]) * scale if r else None
    except (KeyError, ValueError):
        return None


def runs(root, setting, base):
    """[(etiquette, couleur, style, marqueur, [courbes])] dans l'ordre fixe des series."""
    out = []
    if setting in ("ablation", "ablation_k"):
        d = root / ("results_ablation_raw_k" if setting == "ablation_k" else "results_ablation_raw_khalf")
        for mode, lab, color, ls, mk in ABLATION:
            curves = []
            for t in sorted(d.glob(f"{base}_{mode}_s*.trace")):
                curves.append((read_trace(t), end_time(t.with_suffix(".csv"), "time_ms", 1e-3)))
            out.append((lab, color, ls, mk, curves))
        return out
    suffix = "raw_k" if setting == "raw_k" else "raw_khalf"
    for lab, color, ls, mk in METHODS:
        curves = []
        if lab == "Ours":
            for t in sorted((root / f"results_boolmaxls_{suffix}").glob(f"{base}_boolmaxls_s*.trace")):
                curves.append((read_trace(t), end_time(t.with_suffix(".csv"), "time_ms", 1e-3)))
        elif lab == "NuWLS":
            t = root / f"results_nuwls_greedy_{suffix}" / f"{base}_nuwls.trace"
            if t.exists():
                bench = t.with_name(f"{base}_bench.csv")
                curves.append((read_trace(t), end_time(bench, "time_ms", 1e-3)))
        else:
            algo = {"GreConD": "grecond", "GreEss": "greess", "MDL": "mdlgrecond"}[lab]
            t = root / f"results_{algo}_{suffix}" / f"{base}_{algo}_free.trace"
            if t.exists():
                curves.append((read_trace(t), end_time(t.with_suffix(".csv"), "total_time_s", 1.0)))
        out.append((lab, color, ls, mk, curves))
    return out


def plot_one(base, k, series, setting):
    fig, ax = plt.subplots(figsize=(3.4, 2.3))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    drawn, y_max = 0, 0
    for lab, color, ls, mk, curves in series:
        for i, (pts, t_end) in enumerate(curves):
            if not pts:
                continue
            xs = [max(t, T_MIN) for t, _ in pts]
            ys = [e for _, e in pts]
            if t_end and t_end > xs[-1]:
                xs.append(t_end)
                ys.append(ys[-1])
            y_max = max(y_max, max(ys))
            alpha = 1.0 if i == 0 else 0.45
            ax.step(xs, ys, where="post", color=color, linestyle=ls, linewidth=1.4, alpha=alpha,
                    marker=mk, markevery=[len(xs) - 1], markersize=4.5,
                    markerfacecolor=color, markeredgecolor=SURFACE, markeredgewidth=0.8,
                    label=lab if i == 0 else None)
            drawn += 1
    ax.set_xscale("log")
    # 0 visible et ordres de grandeur lisibles pour les erreurs finales faibles.
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(bottom=-0.2, top=max(2.0, y_max * 1.8))
    ax.set_xlabel("time (s)", color=INK2, fontsize=8)
    ax.set_ylabel("error", color=INK2, fontsize=8)
    title = f"{LABELS.get(base, base[0].upper() + base[1:])}, k = {k}"
    ax.set_title(title, color=INK, fontsize=9, loc="left")
    ax.tick_params(colors=INK2, labelsize=7, length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.grid(True, which="major", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    if drawn:
        leg = ax.legend(fontsize=7, frameon=False, labelcolor=INK, loc="best")
        for line in leg.get_lines():
            line.set_linewidth(1.6)
    else:
        ax.text(0.5, 0.5, "no trace", transform=ax.transAxes, ha="center", color=INK2, fontsize=8)
    fig.tight_layout()
    return fig, drawn


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("setting", choices=["raw_k", "raw_khalf", "ablation", "ablation_k"])
    p.add_argument("--root", type=Path, default=here)
    p.add_argument("--out", type=Path, default=here / "figures")
    p.add_argument("--datasets", help="liste separee par des virgules (defaut : toutes)")
    args = p.parse_args()

    lst = "datasets_real.csv" if args.setting in ("raw_k", "ablation_k") else "datasets_real_khalf.csv"
    todo = read_list(args.root / lst)
    if args.datasets:
        keep = set(args.datasets.split(","))
        todo = [(b, k) for b, k in todo if b in keep]
    (args.out / args.setting).mkdir(parents=True, exist_ok=True)
    empty = []
    with PdfPages(args.out / f"{args.setting}.pdf") as pdf:
        for base, k in todo:
            fig, drawn = plot_one(base, k, runs(args.root, args.setting, base), args.setting)
            if not drawn:
                empty.append(base)
            fig.savefig(args.out / args.setting / f"{base}.pdf")
            pdf.savefig(fig)
            plt.close(fig)
    print(f"-> {args.out / args.setting}.pdf ({len(todo)} pages) et {args.out / args.setting}/<matrice>.pdf")
    if empty:
        print(f"[warn] aucune trace pour : {', '.join(empty)}", file=sys.stderr)


if __name__ == "__main__":
    main()
