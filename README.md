# BMF_LNS

Boolean Matrix Factorization par recherche locale pondérée (WLS style
NuWLS) + Large Neighborhood Search MaxSAT.

## Problème

Étant donné `M ∈ {0, 1, -1}^(m×n)` (avec `-1` = valeur manquante ignorée)
et `k`, trouver `A ∈ {0,1}^(m×k)`, `B ∈ {0,1}^(k×n)` minimisant :

```
err(A, B) = #{ (i, j) : M[i,j] ≠ -1 ∧ M[i,j] ≠ (A X B)[i,j] }
```

où `X` est le produit matriciel booléen.

## Méthode : opt4-wls = GI + (WLS + LNS partial)*

| Brique          | Implémentation                        |
|-----------------|---------------------------------------|
| **GI**          | `BMFLocalSearch::initialize_greedy`   |
| **WLS**         | `BMFLocalSearch::solve_weighted`      |
| **LNS partial** | `BMF::lns_step_partial`               |

WLS spécialisée BMF Partial MaxSAT avec les structures du PDF du
superviseur (`nbcover`, `rows_A`, `cols_B`, `pos_score`) :

- flip incrémental `O(|cols_B[l]| · k)`
- **BMS K=10** pour choisir un flip améliorant
- **double flip** à la sortie de minimum local (conditions du PDF)
- **Weighting-PMS** (hard prioritaire) + **force flip**
- `smooth_prob` auto : `0.01` si `n_vars < 10⁴` sinon `0.0003`

LNS partial : sous-problème MaxSAT sur cellules impactées par `(rows, cols)`
avec voisinage 70% top-k / 30% aléatoire, taille adaptative `nh <- 10`,
`+10` sur stagnation, plafond `min(100, min(m, n))`.

Détails algorithmiques : [`docs/ALGORITHME_BMF_LNS.md`](docs/ALGORITHME_BMF_LNS.md).

## Build

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

Sur Narval : `bash cluster/setup_env.sh` charge les modules (`StdEnv/2023`,
`gcc/12.3`, `cmake/3.27.7`, `python/3.11`), compile EvalMaxSAT2022 puis
`rbac_main`.

## Usage

```bash
./rbac_main <file.csv> <k> [timeout_seconds] [--csv-out <path>]
./rbac_main --help
```

- `timeout=0` : pas de limite, Ctrl+C interrompt proprement.
- `--csv-out <path>` : écrit une ligne CSV avec les métriques (schéma en §Sortie).

Exemples :

```bash
./rbac_main ../data/Real/SimplUni/zoo.csv 25 10
./rbac_main ../data/Real/SimplUni/mushroom.csv 51 60 --csv-out /tmp/mushroom.csv
```

## Sortie

Tableau stdout à la fin de l'exécution + CSV optionnel avec 15 colonnes :

```
method, dataset, k, init_errors,
phase1_errors, phase1_last_improve_ms,
phase2_errors, phase2_last_improve_ms,
final_errors, iterations, time_ms,
wls_total_flips, wls_total_time_ms, wls_flips_per_sec,
verif_errors
```

- `phase1` = après greedy + 1re WLS ; `phase2` = fin de la boucle WLS+LNS.
- `wls_*` : débit cumulé de la WLS (utile pour diagnostiquer un goulot).
- `verif_errors` : `-1` non vérifié, `0` `verif.py` d'accord, `>0` désaccord.


## Arborescence

```
src/
  main.cpp                  CLI + pipeline opt4-wls
  BMF.hpp / BMF.cpp         LNS partial (MaxSAT sur voisinage)
  BMFLocalSearch.hpp/.cpp   GI + WLS (structures PDF, flip incremental,
                            double flip, Weighting-PMS)
  SATSolver.hpp / .cpp      Interface EvalMaxSAT2022 (IPAMIR)
  Matrix.hpp                Matrice binaire
  CSVMatrixLoader.hpp       Chargement CSV
  verif.py                  Verification Python (A * B vs M)
data/
  Real/                     Datasets reels (SimplUni, SimplExi, ...)
  Synthetic/                Datasets synthetiques
third_party/
  EvalMaxSAT2022/           Solveur MaxSAT (IPAMIR)
```

## Dépendances

- C++17 (Clang ou GCC ≥ 12)
- CMake ≥ 3.16
- EvalMaxSAT2022 (compilé dans `third_party/`)
- zlib, pthread
- Python 3 + NumPy (pour `verif.py`)
