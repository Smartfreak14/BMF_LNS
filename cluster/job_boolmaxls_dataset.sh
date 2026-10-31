#!/bin/bash
# 1 tache SLURM = 1 dataset x NSEEDS seeds pour BoolMaxLS.
#   SEED_PER_TASK=1 : 1 tache = 1 (dataset, graine) ; la tache t traite la ligne (t-1)/NSEEDS+1
#                     et la graine (t-1)%NSEEDS+1 (tableau de N x NSEEDS taches).
#   INIT_ALGO=greess : BoolMaxLS part des facteurs de cet algorithme (FromBelowBaselines, meme k,
#                     inconnues free) au lieu du glouton ; son temps est retire du budget.
#   BOOLMAXLS_ARGS  : options ajoutees a chaque execution (ex. "--ablation no-double-flip").

#SBATCH --account=def-florenta
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/boolmaxls_%A_%a.out
#SBATCH --error=logs/boolmaxls_%A_%a.err

set -euo pipefail

if [ -z "${PROJECT_ROOT:-}" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_ROOT="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
fi
: "${PROJECT_ROOT:?PROJECT_ROOT non defini}"

CPP_ROOT="$PROJECT_ROOT"
DATA_BASE="${DATA_BASE:-Real}"
RUN_NAME="${RUN_NAME:-boolmaxls}"
NSEEDS="${NSEEDS:-3}"
TIMEOUT_OVERRIDE="${TIMEOUT_OVERRIDE:-}"
REGIME="${REGIME:-raw}"

DATASETS_FILE="${DATASETS_FILE:-$CPP_ROOT/cluster/datasets_real.csv}"
RBAC_BIN="$CPP_ROOT/build/BoolMaxLS"
DATA_DIR="$CPP_ROOT/data/$DATA_BASE"
RESULTS_DIR="$CPP_ROOT/cluster/results_${RUN_NAME}"

mkdir -p "$RESULTS_DIR" "$CPP_ROOT/cluster/logs"

module --force purge
module load StdEnv/2023 gcc/12.3 python/3.11

TASK_ID="${SLURM_ARRAY_TASK_ID:-}"
[ -n "$TASK_ID" ] || exit 1
SEEDS=$(seq 1 "$NSEEDS")
if [ "${SEED_PER_TASK:-0}" = 1 ]; then
    SEEDS=$(( (TASK_ID - 1) % NSEEDS + 1 ))
    TASK_ID=$(( (TASK_ID - 1) / NSEEDS + 1 ))
fi

LINE=$(grep -v '^#' "$DATASETS_FILE" | tail -n +2 | sed -n "${TASK_ID}p")
[ -n "$LINE" ] || exit 1

IFS=',' read -r FILE K TOUT <<< "$LINE"
FILE="$(echo "$FILE" | tr -d '[:space:]')"
K="$(echo "$K" | tr -d '[:space:]')"
TOUT="$(echo "$TOUT" | tr -d '[:space:]')"
[ -n "$TIMEOUT_OVERRIDE" ] && TOUT="$TIMEOUT_OVERRIDE"

CSV_PATH="$DATA_DIR/$FILE"
[ -f "$CSV_PATH" ] || { echo "ERREUR : $CSV_PATH absent"; exit 1; }

BASE="${FILE%.csv}"
read -r -a EXTRA_ARGS <<< "${BOOLMAXLS_ARGS:-}"
echo "=== $FILE k=$K timeout=${TOUT}s graines=$(echo $SEEDS)${INIT_ALGO:+ depart=$INIT_ALGO}${BOOLMAXLS_ARGS:+ options=$BOOLMAXLS_ARGS} ==="

for SEED in $SEEDS; do
    OUT_CSV="$RESULTS_DIR/${BASE}_boolmaxls_s${SEED}.csv"
    TRACE="$RESULTS_DIR/${BASE}_boolmaxls_s${SEED}.trace"
    LOG="$RESULTS_DIR/${BASE}_boolmaxls_s${SEED}.log"
    BUDGET_S="$TOUT"
    INIT_ARGS=()
    if [ -n "${INIT_ALGO:-}" ]; then
        INIT="$RESULTS_DIR/${BASE}_${INIT_ALGO}_init_s${SEED}"
        "$CPP_ROOT/build/FromBelowBaselines" "$INIT_ALGO" "$CSV_PATH" "$K" --missing free \
            --time-limit "$TOUT" --out-dir "$INIT" > "$INIT.log" 2>&1 || true
        [ -f "$INIT/${BASE}_${INIT_ALGO}_free.csv" ] || { echo "ERREUR : $INIT_ALGO sans resultat (voir $INIT.log)"; continue; }
        # total_time_s de l'algorithme, arrondi a la seconde superieure.
        T_INIT=$(awk -F, 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "total_time_s") c = i; next}
                          {print int($c) + ($c > int($c)); exit}' "$INIT/${BASE}_${INIT_ALGO}_free.csv")
        BUDGET_S=$(( TOUT - T_INIT ))
        INIT_ARGS=(--init-factors "$INIT/${BASE}_${INIT_ALGO}_free_k${K}_A.csv" "$INIT/${BASE}_${INIT_ALGO}_free_k${K}_B.csv")
        echo "    s$SEED : $INIT_ALGO en ${T_INIT}s, budget BoolMaxLS ${BUDGET_S}s"
    fi
    "$RBAC_BIN" "$CSV_PATH" "$K" "$BUDGET_S" \
        --seed "$SEED" --regime "$REGIME" ${INIT_ARGS[@]+"${INIT_ARGS[@]}"} ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} \
        --csv-out "$OUT_CSV" --trace-out "$TRACE" > "$LOG" 2>&1 || true
    [ -z "${INIT_ALGO:-}" ] || rm -f "$INIT"/*_A.csv "$INIT"/*_B.csv "$INIT"/*_M.csv
done

echo "Finished : $(date -Iseconds)"
