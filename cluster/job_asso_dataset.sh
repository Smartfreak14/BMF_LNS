#!/bin/bash
# 1 tache SLURM = 1 dataset x NSEEDS seeds pour ASSO (Miettinen 2008).
# Utilise en mode array par submit_asso_real.sh / submit_asso_table2.sh.
#
# Env attendus :
#   PROJECT_ROOT      : racine du projet
#   DATASETS_FILE     : CSV filename,k,timeout
#   DATA_BASE         : "Real" (relatif a data/)
#   RUN_NAME          : suffixe pour results_<RUN_NAME>/
#   NSEEDS            : nombre de graines (defaut 3)
#   TIMEOUT_OVERRIDE  : timeout en secondes (defaut : utilise csv)
#   REGIME            : etiquette 'raw'/'simpluni'/etc

#SBATCH --account=def-florenta
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/asso_%A_%a.out
#SBATCH --error=logs/asso_%A_%a.err

set -euo pipefail

if [ -z "${PROJECT_ROOT:-}" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_ROOT="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
fi
: "${PROJECT_ROOT:?PROJECT_ROOT non defini}"

CPP_ROOT="$PROJECT_ROOT"
DATA_BASE="${DATA_BASE:-Real}"
RUN_NAME="${RUN_NAME:-asso}"
NSEEDS="${NSEEDS:-3}"
TIMEOUT_OVERRIDE="${TIMEOUT_OVERRIDE:-}"
REGIME="${REGIME:-raw}"
TAU_LIST="${TAU_LIST:-0.3,0.5,0.7,0.9}"

DATASETS_FILE="${DATASETS_FILE:-$CPP_ROOT/cluster/datasets_real.csv}"
DATA_DIR="$CPP_ROOT/data/$DATA_BASE"
RESULTS_DIR="$CPP_ROOT/cluster/results_${RUN_NAME}"
LOGS_DIR="$CPP_ROOT/cluster/logs"

mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

module --force purge
module load StdEnv/2023 python/3.11 scipy-stack/2024a
# scipy-stack fournit numpy + scipy + pandas etc.

TASK_ID="${SLURM_ARRAY_TASK_ID:-}"
[ -n "$TASK_ID" ] || { echo "ERREUR : SLURM_ARRAY_TASK_ID absent"; exit 1; }

LINE=$(grep -v '^#' "$DATASETS_FILE" | tail -n +2 | sed -n "${TASK_ID}p")
[ -n "$LINE" ] || { echo "ERREUR : ligne $TASK_ID absente de $DATASETS_FILE"; exit 1; }

IFS=',' read -r FILE K TOUT <<< "$LINE"
FILE="$(echo "$FILE" | tr -d '[:space:]')"
K="$(echo "$K" | tr -d '[:space:]')"
TOUT="$(echo "$TOUT" | tr -d '[:space:]')"
[ -n "$TIMEOUT_OVERRIDE" ] && TOUT="$TIMEOUT_OVERRIDE"

CSV_PATH="$DATA_DIR/$FILE"
[ -f "$CSV_PATH" ] || { echo "ERREUR : $CSV_PATH absent"; exit 1; }

BASE="${FILE%.csv}"

echo "============================================================"
echo " Task $TASK_ID : $FILE  k=$K  timeout=${TOUT}s  nseeds=$NSEEDS"
echo " Regime : $REGIME"
echo " tau_list : $TAU_LIST (budget partage equitablement)"
echo "============================================================"

for SEED in $(seq 1 $NSEEDS); do
    OUT_CSV="$RESULTS_DIR/${BASE}_asso_s${SEED}.csv"
    TRACE="$RESULTS_DIR/${BASE}_asso_s${SEED}.trace"
    LOG="$RESULTS_DIR/${BASE}_asso_s${SEED}.log"
    echo ">>> seed=$SEED"
    python "$CPP_ROOT/scripts/baselines/asso.py" \
        "$CSV_PATH" "$K" "$TOUT" "$OUT_CSV" \
        --trace-out "$TRACE" --seed "$SEED" --regime "$REGIME" \
        --tau-list "$TAU_LIST" \
        > "$LOG" 2>&1 || echo "[warn] seed $SEED failed (voir $LOG)"
done

echo "Finished : $(date -Iseconds)"
