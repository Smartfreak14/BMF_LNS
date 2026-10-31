#!/usr/bin/env python3
"""
Tableau LaTeX au format du tableau principal de l'article (table_main_article.tex), complete
par GreConD, GreEss et MDLGreConD, avec les resultats des campagnes de launch_rerun.sh.

Sources par colonne :
  Ours           : boolmaxls_simpluni / _simplexi / _raw_k / _raw_khalf (mediane des graines) ;
                   a gauche, plus petit rang atteint avec boolmaxls_simpluni_oldk / _simplexi_oldk
  NuWLS (brutes) : nuwls_greedy_raw_k / _raw_khalf (NuWLS-c demarre du glouton), completees par
                   les campagnes _large pour les encodages au-dessus de leur limite
  Sur les brutes, une campagne <nom>_final (configuration finale) remplace <nom> et ses _large.
  Ours compte l'erreur atteinte au budget (--budget, d'apres la trace) : un appel LNS qui deborde
  ne profite pas du temps en plus. Les executions de <campagne>_redo remplacent celles de meme
  nom dans <campagne>.
  GreConD, GreEss, MDL : <algo>_simpluni / _simplexi / _raw_k / _raw_khalf
  autres colonnes (Opti, Kol., ASSO, NuWLS sur les simplifiees) : tableau de l'article
Une colonne dont la campagne n'existe pas encore garde la valeur de l'article ; une matrice
sans resultat dans une campagne existante est marquee "?".

  bloc gauche : rang (temps) de la factorisation exacte, "--" si non atteinte en 30 min
                (MDLGreConD s'arrete par son critere MDL et n'y figure pas)
  bloc droit  : erreur a k et floor(k/2), "--" = budget epuise ou formule non lancee
Gras recalcule sur toutes les colonnes : plus petit rang puis plus petit temps (gauche),
plus petite erreur, egalites toutes en gras (droite).

--diff : en rouge, les valeurs relancees qui different de l'ancien tableau (rang ou erreur ;
les temps ne sont pas compares) ; double dague si la relance a vise un autre rang que lui.
Les listes de k (datasets*.csv) sont lues a cote de ce script, les resultats sous --root.

Usage (depuis cluster/) :
    python make_frombelow_table.py [-o table_results.tex] [--missing free] [--article fichier.tex]
                                   [--root dossier_des_results] [--diff] [--budget 1800]
"""
import argparse, csv, math, re, statistics, sys
from pathlib import Path

NEW = [("grecond", "GreConD"), ("greess", "GreEss"), ("mdlgrecond", "MDL")]
EXACT = [("grecond", "GreConD"), ("greess", "GreEss")]
LABELS = {"americas_small": "Americas\\_small", "chess_krvskp": "Chess",
          "studentPerformance": "StudentPerf.", "ThoraricSurgery": "Thoracic"}
UNIT = {"ms": 1e-3, "s": 1.0, "m": 60.0, "h": 3600.0}
UNKNOWN = ("?", None, "")
DASH = ("--", None, "")
# NuWLS non lance : encodage au-dela de MAX_CLAUSES (launch_rerun.sh), MAX_CLAUSES_LARGE pour _large.
TOO_LARGE = "${}^{\\S}$"
# Bloc gauche : campagnes BoolMaxLS et leurs listes de rangs ; les listes _oldk ajoutent les rangs
# de l'ancien tableau absents des listes principales.
EXACT_RUNS = {"su": [("boolmaxls_simpluni", "datasets.csv"),
                     ("boolmaxls_simpluni_oldk", "datasets_simpluni_oldk.csv")],
              "se": [("boolmaxls_simplexi", "datasets_simplexi.csv"),
                     ("boolmaxls_simplexi_oldk", "datasets_simplexi_oldk.csv")]}


def label(base):
    return LABELS.get(base, base[0].upper() + base[1:])


def fmt_time(t):
    if t < 0.01:
        return f"{t * 1e3:.1f}ms"
    if t < 1:
        return f"{t * 1e3:.0f}ms"
    if t < 10:
        return f"{t:.1f}s"
    if t < 60:
        return f"{t:.0f}s"
    if t < 600:
        return f"{t / 60:.1f}m"
    if t < 3600:
        return f"{t / 60:.0f}m"
    return f"{t / 3600:.1f}h"


def shown(t):
    """Temps tel qu'affiche, pour departager les rangs egaux comme dans l'article."""
    m = re.match(r"^([\d.]+)(ms|s|m|h)$", fmt_time(t))
    return float(m.group(1)) * UNIT[m.group(2)]


def fmt_num(x):
    return str(int(x)) if x == int(x) else f"{x:.1f}"


def sci(n):
    """n arrondi par defaut a un chiffre apres la virgule, en notation scientifique LaTeX."""
    e = len(str(n)) - 1
    return f"{math.floor(n / 10 ** (e - 1)) / 10:.1f} \\cdot 10^{{{e}}}"


def read_k(path):
    rows = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")][1:]
    return {r.split(",")[0].strip()[:-4]: int(r.split(",")[1]) for r in rows}


# Une cellule = (texte sans gras, cle de comparaison ou None, marque LaTeX gardee hors gras).
def frombelow_cells(root, campaign, missing, exact):
    d = root / f"results_{campaign}"
    if not d.is_dir():
        return None
    cells = {}
    # Sans cellule inconnue, seul free est lance (zero et one lui sont identiques).
    for mode in dict.fromkeys(["free", missing]):
        for f in d.glob(f"*_{mode}.csv"):
            r = next(csv.DictReader(f.open()), None)
            if r is None:
                continue
            if r["stop"].startswith("timeout") or (exact and r["error"] != "0"):
                cells[r["dataset"]] = DASH
            elif exact:
                t = float(r["time_s"])
                cells[r["dataset"]] = (f"{r['factors']} ({fmt_time(t)})", (int(r["factors"]), shown(t)), "")
            else:
                cells[r["dataset"]] = (r["error"], (int(r["error"]),), "")
    return cells


def at_budget(f, final, budget_ms):
    """Erreur au budget d'apres la trace : un appel LNS qui finit apres le budget ne compte pas.
    Sans trace (tache tuee avant d'ecrire son tampon), l'erreur du CSV."""
    tr = f.with_suffix(".trace")
    pts = [l.split() for l in tr.read_text().splitlines()[1:]] if tr.exists() else []
    errs = [int(p[1]) for p in pts if len(p) == 2 and float(p[0]) <= budget_ms]
    return errs[-1] if errs else final


def boolmaxls_runs(root, campaign):
    """CSV des executions de results_<campagne> ; ceux de results_<campagne>_redo remplacent
    les executions de meme nom."""
    files = {f.name: f for f in (root / f"results_{campaign}").glob("*_boolmaxls_s*.csv")}
    files.update({f.name: f for f in (root / f"results_{campaign}_redo").glob("*_boolmaxls_s*.csv")})
    return [files[n] for n in sorted(files)]


def boolmaxls_cells(root, campaign, exact, budget_ms):
    d = root / f"results_{campaign}"
    if not d.is_dir():
        return None
    runs = {}
    for f in boolmaxls_runs(root, campaign):
        r = next(csv.DictReader(f.open()), None)
        if r:
            runs.setdefault(f.name.rsplit("_boolmaxls_s", 1)[0], []).append((f, r))
    cells = {}
    for base, rs in runs.items():
        if exact:
            t = [float(r["time_to_zero_ms"]) / 1e3 for _, r in rs
                 if int(r["final_errors"]) == 0 and 0 <= float(r["time_to_zero_ms"]) <= budget_ms]
            if len(t) < len(rs):
                cells[base] = DASH
            else:
                k, tm = int(rs[0][1]["k"]), statistics.median(t)
                cells[base] = (f"{k} ({fmt_time(tm)})", (k, shown(tm)), "")
        else:
            e = statistics.median(at_budget(f, int(r["final_errors"]), budget_ms) for f, r in rs)
            cells[base] = (fmt_num(e), (e,), "")
    return cells


def nuwls_cells(root, campaign, sizes):
    """Mediane des graines (<base>_s<graine>_bench.csv), ou graine unique (<base>_bench.csv).
    sizes : rempli avec le nombre de clauses des encodages non lances."""
    d = root / f"results_{campaign}"
    if not d.is_dir():
        return None
    runs = {}
    for f in d.glob("*_bench.csv"):
        r = next((x for x in csv.DictReader(f.open()) if x["solver"] == "nuwls"), None)
        if r:
            base = re.sub(r"_s\d+$", "", f.name[:-len("_bench.csv")])
            runs.setdefault(base, []).append(r)
    cells = {}
    for base, rs in runs.items():
        costs = [int(r["best_cost"]) for r in rs if r["best_cost"].isdigit()]
        if len(costs) == len(rs):
            e = statistics.median(costs)
            cells[base] = (fmt_num(e), (e,), "")
        elif rs[0].get("note", "").startswith("too_large"):
            cells[base] = ("--", None, TOO_LARGE)
            sizes[base] = int(rs[0]["note"].split("=", 1)[1])
        else:
            cells[base] = DASH
    return cells


def article_cell(text, left):
    text = text.strip()
    mark = ""
    if "$" in text:
        i = text.index("$")
        text, mark = text[:i].strip(), text[i:]
    if text.startswith("\\textbf{") and text.endswith("}"):
        text = text[len("\\textbf{"):-1].strip()
    if left:
        m = re.match(r"^(\d+)\s*\(([\d.]+)\s*(ms|s|m|h)\)$", text)
        key = (int(m.group(1)), float(m.group(2)) * UNIT[m.group(3)]) if m else None
    else:
        key = (int(text),) if text.isdigit() else None
    return (text, key, mark)


def parse_article(path):
    body = path.read_text().split("\\midrule", 1)[1].split("\\bottomrule", 1)[0]
    rows = {}
    for line in body.splitlines():
        line = line.strip()
        if "&" not in line:
            continue
        if line.endswith("\\\\"):
            line = line[:-2]
        c = [x.strip() for x in line.split("&")]
        if len(c) != 16:
            sys.exit(f"ligne de {len(c)} colonnes (16 attendues) : {line[:60]}")
        rows[c[0]] = {"su": [article_cell(x, True) for x in c[1:4]],
                      "se": [article_cell(x, True) for x in c[4:7]],
                      "k": c[7],
                      "rk": [article_cell(x, False) for x in c[8:12]],
                      "rk2": [article_cell(x, False) for x in c[12:16]]}
    return rows


def bold(cells):
    keys = [c[1] for c in cells if c[1] is not None]
    best = min(keys) if keys else None
    return [(f"\\textbf{{{t}}}" if k is not None and k == best else t) + m for t, k, m in cells]


def value(cell):
    """Ce qui est compare a l'ancien tableau : rang (gauche) ou erreur (droite), sinon "--" / "?"."""
    return cell[1][0] if cell[1] is not None else cell[0]


def red(s, on):
    return f"\\textcolor{{red}}{{{s}}}" if on else s


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=here)
    p.add_argument("--missing", default="free", choices=["free", "zero", "one"])
    p.add_argument("--article", type=Path, default=here / "table_main_article.tex")
    p.add_argument("-o", "--output", type=Path)
    p.add_argument("--diff", action="store_true")
    p.add_argument("--budget", type=float, default=1800, help="budget par execution (s)")
    args = p.parse_args()

    root = args.root
    budget_ms = args.budget * 1000
    ks = read_k(here / "datasets_real.csv")
    ks2 = read_k(here / "datasets_real_khalf.csv")
    art = parse_article(args.article)
    m = re.search(r"\\label\{([^}]*)\}", args.article.read_text())
    tab_label = m.group(1) if m else "tab:main-frombelow"
    new = {(a, c): frombelow_cells(root, f"{a}_{c}", args.missing, c.startswith("simpl"))
           for a, _ in NEW for c in ("raw_k", "raw_khalf", "simpluni", "simplexi")}
    warnings = []
    # (bloc, colonne dans le bloc) -> cellules issues des relances, et campagnes lues.
    rerun, used = {}, {}
    # Bloc gauche : plus petit rang atteint, puis plus petit temps ; tried = rangs relances par matrice.
    # Une campagne complementaire n'est lue qu'avec sa campagne principale.
    tried = {}
    for b, runs in EXACT_RUNS.items():
        rerun[(b, 0)], used[(b, 0)], tried[b] = None, [], {}
        for campaign, lst in runs:
            cells = boolmaxls_cells(root, campaign, True, budget_ms)
            if cells is None:
                continue
            if rerun[(b, 0)] is None and campaign != runs[0][0]:
                warnings.append(f"results_{campaign} ignore : results_{runs[0][0]} absent sous {root}")
                continue
            k_run = read_k(here / lst)
            merged = rerun[(b, 0)] = rerun[(b, 0)] or {}
            used[(b, 0)].append(campaign)
            for base, c in cells.items():
                tried[b].setdefault(base, []).append(k_run.get(base))
                cur = merged.get(base)
                if cur is None or (c[1] is not None and (cur[1] is None or c[1] < cur[1])):
                    merged[base] = c
    # Brutes : la configuration finale (_final, plusieurs graines) remplace la premiere relance.
    for b, campaign in (("rk", "boolmaxls_raw_k"), ("rk2", "boolmaxls_raw_khalf")):
        if (root / f"results_{campaign}_final").is_dir():
            campaign += "_final"
        rerun[(b, 0)], used[(b, 0)] = boolmaxls_cells(root, campaign, False, budget_ms), [campaign]
        if (root / f"results_{campaign}_redo").is_dir():
            used[(b, 0)].append(campaign + "_redo")
    # NuWLS : une campagne _large remplace les encodages trop gros pour la limite de la premiere.
    sizes = {}
    for b, campaign in (("rk", "nuwls_greedy_raw_k"), ("rk2", "nuwls_greedy_raw_khalf")):
        sizes[b] = {}
        if (root / f"results_{campaign}_final").is_dir():
            rerun[(b, 1)] = nuwls_cells(root, campaign + "_final", sizes[b])
            used[(b, 1)] = [campaign + "_final"]
            continue
        cells = nuwls_cells(root, campaign, sizes[b])
        large = nuwls_cells(root, campaign + "_large", sizes[b])
        used[(b, 1)] = [campaign]
        if large is not None and cells is None:
            warnings.append(f"results_{campaign}_large ignore : results_{campaign} absent sous {root}")
        elif large is not None:
            cells.update(large)
            used[(b, 1)].append(campaign + "_large")
        rerun[(b, 1)] = cells

    nl, nr = 3 + len(EXACT), 4 + len(NEW)
    head_left = " & ".join(["Ours", "NuWLS", "Opti"] + [n for _, n in EXACT])
    head_right = " & ".join(["Ours", "NuWLS", "Kol.", "ASSO"] + [n for _, n in NEW])
    # Colonnes : matrice, Simpli-forall, Simpli-exists0, k, brutes a k, floor(k/2), brutes a floor(k/2).
    c1, c2, c3 = 1 + nl, 1 + 2 * nl, 3 + 2 * nl
    lines = [
        "\\begin{table*}[h]", "\\centering",
        "\\setlength{\\tabcolsep}{0.2pt}\\renewcommand{\\arraystretch}{0.8}", "\\scriptsize",
        "\\resizebox{\\textwidth}{!}{%",
        f"\\begin{{tabular}}{{l {'r' * nl} {'r' * nl} r {'r' * nr} r {'r' * nr}}}", "\\toprule",
        f" & \\multicolumn{{{nl}}}{{c}}{{$\\text{{Simpli}}^{{\\forall}}$: $k$ (time)}} "
        f"& \\multicolumn{{{nl}}}{{c}}{{$\\text{{Simpli}}^{{\\exists}}_{{0}}$: $k$ (time)}} & "
        f"& \\multicolumn{{{nr}}}{{c}}{{raw, rank $k$}} & "
        f"& \\multicolumn{{{nr}}}{{c}}{{raw, rank $\\lfloor k/2 \\rfloor$}} \\\\",
        f"\\cmidrule(lr){{2-{c1}}}\\cmidrule(lr){{{c1 + 1}-{c2}}}"
        f"\\cmidrule(lr){{{c3}-{c3 + nr - 1}}}\\cmidrule(lr){{{c3 + nr + 1}-{c3 + 2 * nr}}}",
        f"\\textsc{{Dataset}} & {head_left} & {head_left} & $k$ & {head_right} "
        f"& $\\lfloor k/2 \\rfloor$ & {head_right} \\\\",
        "\\midrule"]

    COLS = {"su": "Simpli-forall", "se": "Simpli-exists0", "rk": "brutes k", "rk2": "brutes k/2"}
    diffs, dagger = [], False
    for base in sorted(ks, key=str.lower):
        name = label(base)
        old = art.get(name)
        if old is None:
            sys.exit(f"{name} absent du tableau de l'article {args.article}")
        blocks = {b: list(old[b]) for b in ("su", "se", "rk", "rk2")}
        changed, notes = set(), []
        for (b, i), cells in rerun.items():
            if cells is not None:
                c = cells.get(base, UNKNOWN)
                if args.diff and value(c) != value(old[b][i]):
                    changed.add((b, i))
                    diffs.append(f"{name:15s} {COLS[b]:15s} {'Ours' if i == 0 else 'NuWLS':5s} : "
                                 f"{old[b][i][0]} -> {c[0]}{' (non lance)' if c[2] == TOO_LARGE else ''}")
                    k_old = old[b][i][1][0] if b in tried and old[b][i][1] is not None else None
                    if k_old is not None and k_old not in tried[b].get(base, []):
                        c, dagger = (c[0], c[1], c[2] + "${}^{\\ddagger}$"), True
                        notes.append(f"{b} : rang relance {'/'.join(map(str, tried[b].get(base, [])))}, "
                                     f"ancien tableau {k_old}")
                blocks[b][i] = c
        for b, c, algos in (("su", "simpluni", EXACT), ("se", "simplexi", EXACT),
                            ("rk", "raw_k", NEW), ("rk2", "raw_khalf", NEW)):
            for a, _ in algos:
                cells = new[(a, c)]
                blocks[b].append(UNKNOWN if cells is None else cells.get(base, UNKNOWN))
        if old["k"] != str(ks[base]):
            notes.append(f"k article = {old['k']}, k des campagnes = {ks[base]}")
            warnings.append(f"{name} : k article = {old['k']}, datasets_real.csv = {ks[base]}")
        if ks2.get(base) != ks[base] // 2:
            warnings.append(f"{name} : datasets_real_khalf.csv = {ks2.get(base)}, floor(k/2) = {ks[base] // 2}")
        row = {b: [red(s, (b, i) in changed) for i, s in enumerate(bold(blocks[b]))] for b in blocks}
        k_cell = red(str(ks[base]), args.diff and old["k"] != str(ks[base]))
        cells = row["su"] + row["se"] + [k_cell] + row["rk"] + [str(ks2.get(base, "?"))] + row["rk2"]
        note = f"  % {'; '.join(notes)}" if notes else ""
        lines.append(f"{name} & " + " & ".join(cells) + " \\\\" + note)

    reran = [("Ours", any(rerun[(b, 0)] is not None for b in ("su", "se", "rk", "rk2"))),
             ("NuWLS on the raw matrices (NuWLS-c started from the greedy solution)",
              any(rerun[(b, 1)] is not None for b in ("rk", "rk2")))]
    reran = [n for n, ok in reran if ok]
    unknown_txt = ("Unknown cells may be covered and are not counted in the error, as for the other methods. "
                   if args.missing == "free" else
                   f"Unknown cells set to {'0' if args.missing == 'zero' else '1'} for these three methods. ")
    seeds = {f.stem.rsplit("_s", 1)[1] for f in (root / f"results_{used[('rk', 0)][0]}").glob("*_boolmaxls_s*.csv")}
    seeds_txt = (f"median over {len(seeds)} seeds for the randomized methods" if len(seeds) > 1
                 else "one run for the randomized methods")
    diff_txt = ("\\textcolor{red}{Red}: differs from the previous version of this table (rank or error; "
                "times are not compared)"
                + ("; ${}^{\\ddagger}$ re-run at another rank than in that version" if dagger else "") + ". "
                if args.diff else "")
    skipped = [sizes[b][base] for b in ("rk", "rk2") if rerun[(b, 1)]
               for base, c in rerun[(b, 1)].items() if c[2] == TOO_LARGE]
    too_large_txt = (f"{TOO_LARGE} NuWLS not run: Partial MaxSAT encoding of ${sci(min(skipped))}$ clauses or more. "
                     if skipped else "")
    lines += ["\\bottomrule", "\\end{tabular}}",
              "\\caption{Exact factorization (left) and reconstruction error (right) on the 31 benchmarks. "
              "\\emph{Ours} = \\ourmethod{}, \\emph{NuWLS} = NuWLS-c, \\emph{Opti} = OptiBlock, "
              "\\emph{Kol.} = \\citet{kolomvakis2025algorithms}, GreConD = Belohlavek and Vychodil (2010), "
              "GreEss = Belohlavek and Trnecka (2015), MDL = MDLGreConD, Makhalova and Trnecka (2021); "
              "the last three are deterministic (one run, 30-min limit). "
              + unknown_txt
              + (f"Re-run: {' and '.join(reran)}; other values from the previous version of this table. "
                 if reran else "")
              + diff_txt
              + "\\textbf{Left block:} rank and time to reach an exact factorization on the simplified matrices. "
                "``--'' = no zero error within the 30-min limit. ${}^{\\dagger}$ OptiBlock best (not always zero) "
                "under 3h. MDLGreConD stops on its MDL criterion and is not shown. "
                "Bold: smallest rank, ties broken by shortest time. "
                "\\textbf{Right block:} error on the raw matrices at rank $k$ and $\\lfloor k/2 \\rfloor$ "
                f"({seeds_txt}); ``--'' = timeout (30 min). " + too_large_txt +
                "Bold: best per setting, ties all bolded.}",
              f"\\label{{{tab_label}}}", "\\end{table*}"]

    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
        print(f"-> {args.output}", file=sys.stderr)
    else:
        print(text)
    for (b, i), cells in rerun.items():
        col = COLS[b]
        who = "Ours" if i == 0 else "NuWLS"
        src = ("valeurs de l'article" if cells is None else
               f"relance ({len(cells)}/{len(ks)} matrices : {' + '.join(used.get((b, i), []))})")
        print(f"[source] {who:5s} {col:15s} : {src}", file=sys.stderr)
    for (a, c), cells in new.items():
        if c.startswith("simpl") and a == "mdlgrecond":
            continue
        n = 0 if cells is None else len(cells)
        if n < len(ks):
            print(f"[warn] results_{a}_{c} : {n}/{len(ks)} matrices (\"?\" dans le tableau)", file=sys.stderr)
    for w in warnings:
        print(f"[warn] {w}", file=sys.stderr)
    for d in diffs:
        print(f"[diff] {d}", file=sys.stderr)


if __name__ == "__main__":
    main()
