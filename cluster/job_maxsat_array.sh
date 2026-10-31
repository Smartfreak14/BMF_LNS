#!/bin/bash
# 1 tache SLURM = 1 dataset : dump WCNF puis lance les solveurs MaxSAT
# anytime (EvalMaxSAT, NuWLS, SATLike-c, loandra) sur cette formule.
#
# Variables d'environnement (set par submit_maxsat_*.sh) :
#   PROJECT_ROOT   : racine du projet
#   DATASETS_FILE  : CSV filename,k,timeout
#   DATA_BASE      : ex. "Real/SimplUni"
#   RUN_NAME       : suffixe pour results_maxsat_<RUN_NAME>/
#   MAXSAT_TIMEOUT : timeout par solveur (secondes). Defaut 60.
#   SEED_PER_TASK  : 1 = 1 tache par (dataset, graine) parmi NSEEDS (voir job_boolmaxls_dataset.sh),
#                    sorties <base>_s<graine>_* ; sinon une seule graine (SEED).

#SBATCH --account=def-florenta
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/maxsat_%A_%a.out
#SBATCH --error=logs/maxsat_%A_%a.err
#
# --mem=32G : les baselines MaxSAT sur formules BMF grosses (40M+ clauses)
# demandent plusieurs Go. 8G ne suffit pas pour wine, mushroom, etc.

set -euo pipefail

if [ -z "${PROJECT_ROOT:-}" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_ROOT="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
fi
: "${PROJECT_ROOT:?PROJECT_ROOT non defini}"

if [ -f "$PROJECT_ROOT/CMakeLists.txt" ]; then
    CPP_ROOT="$PROJECT_ROOT"
elif [ -f "$PROJECT_ROOT/BoolMaxLStenance/cpp/CMakeLists.txt" ]; then
    CPP_ROOT="$PROJECT_ROOT/BoolMaxLStenance/cpp"
else
    echo "ERREUR : CMakeLists.txt introuvable sous $PROJECT_ROOT"; exit 1
fi

DATA_BASE="${DATA_BASE:-Real/SimplUni}"
RUN_NAME="${RUN_NAME:-maxsat_baseline}"
MAXSAT_TIMEOUT="${MAXSAT_TIMEOUT:-60}"

DATASETS_FILE="${DATASETS_FILE:-$CPP_ROOT/cluster/datasets.csv}"
RBAC_BIN="$CPP_ROOT/build/BoolMaxLS"
DATA_DIR="$CPP_ROOT/data/$DATA_BASE"
RESULTS_DIR="$CPP_ROOT/cluster/results_${RUN_NAME}"
LOGS_DIR="$CPP_ROOT/cluster/logs"

mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

module --force purge
module load StdEnv/2023
module load gcc/12.3
module load python/3.11

TASK_ID="${SLURM_ARRAY_TASK_ID:-}"
if [ -z "$TASK_ID" ] || [ "$TASK_ID" = "4294967294" ]; then
    echo "ERREUR : SLURM_ARRAY_TASK_ID invalide ('$TASK_ID')."; exit 1
fi
SUFFIX=""
if [ "${SEED_PER_TASK:-0}" = 1 ]; then
    : "${NSEEDS:?NSEEDS non defini}"
    export SEED=$(( (TASK_ID - 1) % NSEEDS + 1 ))
    TASK_ID=$(( (TASK_ID - 1) / NSEEDS + 1 ))
    SUFFIX="_s$SEED"
fi

LINE=$(grep -v '^#' "$DATASETS_FILE" | tail -n +2 | sed -n "${TASK_ID}p")
[ -n "$LINE" ] || { echo "ERREUR : pas de ligne $TASK_ID dans $DATASETS_FILE"; exit 1; }

IFS=',' read -r FILE K _IGNORED <<< "$LINE"
FILE="$(echo "$FILE" | tr -d '[:space:]')"
K="$(echo "$K" | tr -d '[:space:]')"

CSV_PATH="$DATA_DIR/$FILE"
[ -f "$CSV_PATH" ] || { echo "ERREUR : dataset $CSV_PATH absent"; exit 1; }

BASE="${FILE%.csv}"
WCNF="$RESULTS_DIR/${BASE}${SUFFIX}.wcnf"
BENCH_CSV="$RESULTS_DIR/${BASE}${SUFFIX}_bench.csv"
OUT_LOG="$RESULTS_DIR/${BASE}${SUFFIX}.log"

echo "============================================================"
echo " Task $TASK_ID  (run=$RUN_NAME)"
echo " Dataset  : $DATA_BASE/$FILE (k=$K, graine ${SEED:-42})"
echo " Timeout  : ${MAXSAT_TIMEOUT}s par solveur"
echo " WCNF     : $WCNF"
echo " Bench    : $BENCH_CSV"
echo " Started  : $(date -Iseconds)"
echo "============================================================"

MAX_CLAUSES="${MAX_CLAUSES:-20000000}"
# NUWLS_INIT=greedy : NuWLS demarre de l'affectation gloutonne de BoolMaxLS.
NUWLS_INIT="${NUWLS_INIT:-}"
INIT_FILE="$RESULTS_DIR/${BASE}${SUFFIX}.init"
INIT_ARGS=()
[ "$NUWLS_INIT" = "greedy" ] && INIT_ARGS=(--init-out "$INIT_FILE")

{
    # Au-dela de MAX_CLAUSES, seul l'en-tete est ecrit (pas de formule de plusieurs Go).
    "$RBAC_BIN" --dump-wcnf "$CSV_PATH" "$K" "$WCNF" --max-clauses "$MAX_CLAUSES" ${INIT_ARGS[@]+"${INIT_ARGS[@]}"}

    NCLAUSES=$(head -5 "$WCNF" | grep "^p wcnf" | awk '{print $4}')
    if [ -n "$NCLAUSES" ] && [ "$NCLAUSES" -gt "$MAX_CLAUSES" ]; then
        echo "[skip] $FILE : nclauses=$NCLAUSES > $MAX_CLAUSES, formule trop grosse"
        echo "solver,best_cost,time_ms,timeout_s,solved,note,first_cost,first_time_s" > "$BENCH_CSV"
        echo "nuwls,N/A,0,$MAXSAT_TIMEOUT,0,too_large_nclauses=$NCLAUSES,," >> "$BENCH_CSV"
    elif [ "$NUWLS_INIT" = "greedy" ]; then
        NUWLS_INIT_FILE="$INIT_FILE" bash "$CPP_ROOT/scripts/bench_maxsat.sh" "$WCNF" "$MAXSAT_TIMEOUT" "$BENCH_CSV"
    else
        bash "$CPP_ROOT/scripts/bench_maxsat.sh" "$WCNF" "$MAXSAT_TIMEOUT" "$BENCH_CSV"
    fi
} > "$OUT_LOG" 2>&1
RC=$?

# Nettoyer le WCNF et l'affectation (peuvent peser plusieurs Go).
rm -f "$WCNF" "$INIT_FILE"

echo "Finished : $(date -Iseconds)  exit=$RC"
exit $RC
