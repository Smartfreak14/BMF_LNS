#!/bin/bash
# 1 tache SLURM = 1 dataset x 5 configurations d'ablation x NSEEDS seeds.
#
# Cinq configurations (chaine cumulative de contributions) :
#   greedy         : GreedyInit seul (benefice de l'init gloutonne)
#   no-weighting   : GI + WLS non pondere, sans DF, sans LNS (structure + flips)
#   wls-no-double-no-lns : GI + WLS pondere, sans DF, sans LNS (isole ponderation)
#   no-lns         : GI + WLS pondere + DF, sans LNS (isole double flip)
#   full           : pipeline complet (isole LNS)
#
# Env : PROJECT_ROOT, DATASETS_FILE, DATA_BASE, RUN_NAME, NSEEDS, REGIME.
#
# Timeout par config = colonne 3 du CSV. Timeout uniforme recommande.
# Budget total = 5 configs x NSEEDS x timeout.
# Avec NSEEDS=3 et timeout=1800s : 5*3*1800 = 27000s = 7.5h (marge OK).

#SBATCH --account=def-florenta
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/ablation_%A_%a.out
#SBATCH --error=logs/ablation_%A_%a.err

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
    echo "ERREUR : CMakeLists.txt introuvable"; exit 1
fi

# DATA_BASE defaut = Real/SimplUni (dataset simplifie de Loom/AV24).
DATA_BASE="${DATA_BASE:-Real/SimplUni}"
RUN_NAME="${RUN_NAME:-ablation}"
NSEEDS="${NSEEDS:-3}"
REGIME="${REGIME:-simpluni}"

DATASETS_FILE="${DATASETS_FILE:-$CPP_ROOT/cluster/datasets_ablation.csv}"
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
[ -n "$TASK_ID" ] || { echo "ERREUR : SLURM_ARRAY_TASK_ID absent"; exit 1; }

LINE=$(grep -v '^#' "$DATASETS_FILE" | tail -n +2 | sed -n "${TASK_ID}p")
[ -n "$LINE" ] || { echo "ERREUR : ligne $TASK_ID absente"; exit 1; }

IFS=',' read -r FILE K TOUT <<< "$LINE"
FILE="$(echo "$FILE" | tr -d '[:space:]')"
K="$(echo "$K" | tr -d '[:space:]')"
TOUT="$(echo "$TOUT" | tr -d '[:space:]')"
[ -n "${TIMEOUT_OVERRIDE:-}" ] && TOUT="$TIMEOUT_OVERRIDE"

CSV_PATH="$DATA_DIR/$FILE"
[ -f "$CSV_PATH" ] || { echo "ERREUR : $CSV_PATH absent"; exit 1; }

BASE="${FILE%.csv}"

echo "============================================================"
echo " Task $TASK_ID : $FILE  k=$K  timeout=${TOUT}s / config  regime=$REGIME"
echo " Data dir : $DATA_DIR"
echo " Results  : $RESULTS_DIR"
echo "============================================================"

# 5 configurations x NSEEDS graines
# Ordre logique : de la plus simple (greedy) a la plus complete (full)
MODES=(greedy no-weighting wls-no-double-no-lns no-lns full)

for MODE in "${MODES[@]}"; do
    for SEED in $(seq 1 $NSEEDS); do
        OUT_CSV="$RESULTS_DIR/${BASE}_${MODE}_s${SEED}.csv"
        TRACE="$RESULTS_DIR/${BASE}_${MODE}_s${SEED}.trace"
        OUT_LOG="$RESULTS_DIR/${BASE}_${MODE}_s${SEED}.log"
        echo ">>> mode=$MODE seed=$SEED"
        "$RBAC_BIN" "$CSV_PATH" "$K" "$TOUT" \
            --ablation "$MODE" --seed "$SEED" \
            --regime "$REGIME" \
            --csv-out "$OUT_CSV" --trace-out "$TRACE" > "$OUT_LOG" 2>&1 || echo "[warn] mode=$MODE seed=$SEED failed"
    done
done

echo "Finished : $(date -Iseconds)"
