# BoolMaxLS

Anytime from-below Boolean Matrix Factorization via MaxSAT-style local search.

Given a binary matrix `M` and a target rank `k`, BoolMaxLS finds two binary factors `A`
(`m × k`) and `B` (`k × n`) whose Boolean product `A x B` covers no 0-entry of `M`
(from-below constraint) and covers as many 1-entries as possible. It combines a
greedy from-below initialisation, a weighted local search on the factor matrices,
a structural double-flip move, and a large-neighborhood search that solves
restricted subproblems exactly with a complete Partial MaxSAT solver.



## Build

Requirements: a C++17 compiler, CMake ≥ 3.16, zlib, pthread, Python 3 (for the
verification script and the ASSO baseline). Gurobi and MATLAB are only needed
to run the Kolomvakis baseline.

```bash
./build.sh
```

This produces `build/BoolMaxLS`.

`third_party/NuWLS` (NuWLS-c) and `third_party/BooleanMF` (Kolomvakis et al.) are git
submodules. After `git submodule update --init`, apply the local changes with
`git -C third_party/NuWLS apply ../patches/NuWLS.patch` (option `-init FILE`, an initial
assignment for the first try, used by the `nuwls_greedy_*` campaigns; no `-static`) and
`git -C third_party/BooleanMF apply ../patches/BooleanMF.patch` (adds `utils/vec.m`).

## Run

```
./build/BoolMaxLS <matrix.csv> <k> [timeout_seconds] [OPTIONS]
```

Input format: comma-separated matrix, `0` / `1` for known cells and an empty
field for undefined (don't-care) entries.

Common options:

```
--ablation MODE     full | greedy | no-weighting | wls-no-double-no-lns |
                    no-lns | no-double-flip
--seed N            RNG seed (default 42)
--csv-out PATH      write a one-line CSV summary (unified schema)
--trace-out PATH    write "time_ms err" per improvement
--regime LABEL      raw | simpluni | simplexi | synth (stored in the CSV)
--init-factors A B  start from these factors (CSV, at most k, e.g. GreEss) instead of the greedy
--lns-time-limit S  solve each exact LNS subproblem in a child process, abandoned after S seconds
                    and at the latest at the end of the budget (inf: only at the end of the
                    budget); without this option, an LNS call has no time limit
--dump-wcnf FILE    dump the Partial MaxSAT encoding of the full problem
--help              show all options
```

Example:

```bash
./build/BoolMaxLS data/Real/zoo.csv 25 60 --ablation full --seed 1 \
    --csv-out zoo_full_s1.csv --trace-out zoo_full_s1.trace
```

Output columns of the CSV: `method, dataset, k, init_errors, phase1_errors,
phase1_last_improve_ms, phase2_errors, phase2_last_improve_ms, final_errors,
iterations, time_ms, wls_total_flips, wls_total_time_ms, wls_flips_per_sec,
verif_errors, m, n, regime, tool, seed, time_to_zero_ms, trace_file`.

`verif_errors` is filled by an independent Python script (`src/verif.py`) that
recomputes `A x B` and reports both the number of covered 0-entries (hard
violations) and the number of uncovered 1-entries. A warning is printed to
`stderr` whenever the returned factorisation is not from-below.


## Baselines

Four external methods can be compared:

- **BoolMaxLS** (this work) — `./build/BoolMaxLS ...`.
- **NuWLS-c** — anytime Partial MaxSAT solver run on the WCNF produced by
  `./build/BoolMaxLS --dump-wcnf`.
- **Kolomvakis et al. 2025** — MATLAB + Gurobi, driver in
  `scripts/baselines/run_kolomvakis.m`.
- **ASSO** 
  `scripts/baselines/asso.py`
  `--tau-list` for the confidence-threshold
  sweep documented in the paper.
- **GreConD, GreEss, MDLGreConD** (deterministic from-below heuristics) —
  `./build/FromBelowBaselines <grecond|greess|mdlgrecond> <matrix.csv> <k>
  [--cuts k,k2] [--missing free|zero|one] [--time-limit s] [--out-dir dir]`.
  Factors are added one at a time, so one run up to `k` also gives the first
  `k2` factors. Unknown cells: `free` (default) lets factors cover them without
  having to, as for the other methods; `zero` / `one` replace them before the
  run. Errors are counted on known cells only.
  `cluster/make_frombelow_table.py` builds the LaTeX results table from the
  `<algo>_raw_k`, `_raw_khalf`, `_simpluni` and `_simplexi` campaigns, the
  BoolMaxLS and NuWLS re-runs, and `cluster/table_main_article.tex` for the
  methods that were not re-run. `--root` points to the directory holding the
  `results_*` folders; `--diff` prints in red the re-run values that differ
  from `table_main_article.tex` and lists them on stderr. The
  `boolmaxls_simpl*_oldk` campaigns add the ranks of that table missing from
  `datasets.csv` / `datasets_simplexi.csv`, and the left block keeps the
  smallest rank reached; the `nuwls_greedy_*_large` campaigns run NuWLS on the
  encodings between `MAX_CLAUSES` and `MAX_CLAUSES_LARGE` (`launch_rerun.sh`).
  On the raw matrices, a `<campaign>_final` run (final configuration, several
  seeds) replaces `<campaign>` and its `_large` complement, and the error of
  BoolMaxLS is read from its trace at the budget (`--budget`), so an LNS call
  that ends after the budget does not count. `cluster/make_synthetic_table.py` builds the table of
  the synthetic instances (time to an exact factorization at the planted rank, from
  `boolmaxls_synthetic`, `nuwls_greedy_synthetic` and `<algo>_synthetic`);
  `cluster/make_greess_start_table.py` compares BoolMaxLS started from the greedy solution and from
  the GreEss factors (`boolmaxls_greess_raw_k*`) with GreEss itself; `cluster/make_ablation_table.py`
  builds the ablation table (GI, no-w, +wt, +DF, full) on the raw matrices from `ablation_raw_k` and
  `ablation_raw_khalf`; `cluster/make_nodf_table.py` compares the full method with the full method
  without the double flip (`boolmaxls_nodf_raw_k`, `_raw_khalf`, same settings as the `_final`
  campaigns) and prints on how many matrices removing it worsens or improves the median error. The
  runs of `results_<campaign>_redo` (e.g. `boolmaxls_raw_k_final_redo`, re-run with
  `--lns-time-limit inf` because SLURM killed the original runs inside an LNS call) replace the runs
  of the same name in `results_<campaign>`. `cluster/status_rerun.sh` reports the progress
  of every campaign and builds `cluster/table_results.tex` once all are done.
  Each submitted array is recorded in `cluster/submitted_jobs.log`, so a
  campaign without logs is reported as not submitted, or as submitted with its
  `sacct` state.
  `cluster/make_dataset_table.py` prints the dataset table (m, n, density,
  known 1s, unknown cells, k, min(m,n)); `cluster/plot_traces.py
  raw_k|raw_khalf|ablation|ablation_k` draws the error-time curves from the `.trace` files
  (matplotlib; `module load scipy-stack` on Narval).

## Cluster campaigns (SLURM)

The `cluster/` directory contains the launch scripts used for the SLURM
campaigns on Compute Canada. All scripts share the helper `_submit_common.sh`
and read a dataset list from `datasets*.csv`.

### Deployment

`cluster/deploy_narval.sh`, run on the local machine, creates the project
directory on the cluster and copies the sources, the `third_party` code, the
scripts, the dataset lists and the data in one `rsync`, never local binaries
or object files, together with a SHA-256 manifest (`SOURCES.sha256`). On the
cluster login node, `cluster/setup_env.sh` builds the project that contains it
(an exported `PROJECT_ROOT` is ignored), stops if the manifest is missing or
does not match, builds EvalMaxSAT2022 if the library is missing or not an ELF
object, rebuilds BoolMaxLS, FromBelowBaselines and NuWLS-c from scratch, and
runs short sanity tests. `cluster/launch_rerun.sh` sets every variable its jobs
read (project root, data directory, NuWLS binary, unknown-cell treatments), so
variables exported in the shell cannot change the campaigns, and it skips a
campaign whose job is still queued. The final-configuration campaigns
(`*_final`, `boolmaxls_synthetic`, `nuwls_greedy_synthetic`,
`boolmaxls_greess_raw_k*`) run one SLURM task per (matrix, seed), with
`NSEEDS_FINAL` seeds (default 3); `boolmaxls_greess_*` starts BoolMaxLS from
the GreEss factors computed in the same task, whose time is taken from the budget.
`boolmaxls_nodf_raw_k*` (`--ablation no-double-flip`) and `boolmaxls_raw_k*_final_redo`
(matrices of `datasets_real*_redo.csv`) work the same way, with `--lns-time-limit inf`: the
error at 30 minutes is the same as without a limit, but every run ends and is verified. The
job passes the options in `BOOLMAXLS_ARGS`, which `launch_rerun.sh` sets for every campaign.

### Available launchers

| Method     | Submit script                                            | Backing job script                |
|------------|----------------------------------------------------------|-----------------------------------|
| BoolMaxLS  | `submit_boolmaxls_real.sh`, `submit_boolmaxls_table2.sh` | `job_boolmaxls_dataset.sh`        |
| ASSO       | `submit_asso_real.sh`, `submit_asso_sweep.sh`, `submit_asso_table2.sh` | `job_asso_dataset.sh` |
| Kolomvakis | `submit_kolomvakis.sh`, `submit_kolomvakis_batch.sh`     | `job_kolomvakis_dataset.sh`       |
| NuWLS-c    | `submit_maxsat_{real,simpluni,simplexi,synthetic,table2}.sh` | `job_maxsat_array.sh`         |
| Ablation   | `submit_ablation.sh`                                     | `job_ablation_array.sh`           |
| GreConD, GreEss, MDLGreConD | `launch_rerun.sh` (`<algo>_raw_k`, `_raw_khalf`, `_simpluni`, `_simplexi`, `_synthetic`) | `job_frombelow_dataset.sh` |

`launch_table2_all.sh` is a master script that runs all four methods on the
Table 2 grid (`--method M --budget T --rank R`).

### Typical invocation

```bash
# One-off (from cluster/):
NSEEDS=11 TIMEOUT_OVERRIDE=1800 bash submit_boolmaxls_real.sh

# Common overrides:
#   NSEEDS            : seeds per (method, dataset)  (default varies)
#   TIMEOUT_OVERRIDE  : per-run wall-clock in seconds
#   TIME_LIMIT        : SLURM --time (HH:MM:SS)
#   MEM_PER_TASK      : SLURM --mem
#   RUN_NAME          : suffix for results_<RUN_NAME>/ and log filenames
```

### Aggregation

Each SLURM run produces one CSV per (dataset, seed) in `cluster/results_<RUN_NAME>/`.
The aggregators consolidate them into a single wide-format CSV per campaign:

```bash
python cluster/aggregate_table2.py    results_<RUN_NAME>/ -o table2_<name>.csv
python cluster/aggregate_ablation.py  results_<RUN_NAME>/ -o ablation_<name>.csv
python cluster/aggregate_kolomvakis.py results_<RUN_NAME>/ -o kol_<name>.csv
python cluster/aggregate_maxsat.py    results_<RUN_NAME>/ -o maxsat_<name>.csv
```

## Data

Benchmarks live under `data/`:

- `data/Real/` — 31 UCI-style matrices (raw).
- `data/Real/SimplUni/` — same matrices preprocessed by the universal
  delegation-relegation simplification of Avellaneda 2024.
- `data/Real/SimplExi/` — the aggressive existential variant.
- `data/Synthetic/` — synthetic matrices `A* x B*` generated per the protocol
  of Avellaneda 2024.

## Verification

For every returned factorisation, `src/verif.py` recomputes `A x B` and prints:

```
Covered 0-entries (hard violations) = 0
Uncovered 1-entries (soft violations, err_M) = 0
From-below = yes
Reconstruction error = 0
```

The C++ driver invokes it after every run (not only on exact success) and
prints a warning to `stderr` if `Covered 0-entries > 0`.
