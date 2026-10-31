#!/bin/bash
# Lance BoolMaxLS en PARALLELE sur 31 datasets Real avec K/2 (datasets_table2.csv).
# Meme structure que submit_boolmaxls_real.sh.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

DATASETS_FILE="${DATASETS_FILE:-$HERE/datasets_table2.csv}"
DATA_BASE="Real"
NSEEDS="${NSEEDS:-3}"
TIMEOUT_OVERRIDE="${TIMEOUT_OVERRIDE:-}"
REGIME="${REGIME:-raw}"
RUN_NAME="${RUN_NAME:-boolmaxls_table2}"

ACCOUNT="${ACCOUNT:-def-florenta}"
MEM="${MEM:-16G}"
TIME_LIMIT="${TIME_LIMIT:-12:00:00}"
ARRAY_LIMIT="${ARRAY_LIMIT:-}"

if [ -d "$HERE/../src" ] && [ -d "$HERE/../data" ]; then
    CPP_ROOT="$(cd "$HERE/.." && pwd)"
else
    echo "ERREUR : structure introuvable"; exit 1
fi

[ -x "$CPP_ROOT/build/BoolMaxLS" ] || { echo "ERREUR : BoolMaxLS absent"; exit 1; }
[ -f "$DATASETS_FILE" ] || { echo "ERREUR : $DATASETS_FILE absent"; exit 1; }
[ -f "$HERE/job_boolmaxls_dataset.sh" ] || { echo "ERREUR : job_boolmaxls_dataset.sh absent"; exit 1; }

n=$(awk 'NR==1 {next} /^#/ {next} NF {x++} END {print x+0}' "$DATASETS_FILE")

mkdir -p "$CPP_ROOT/cluster/logs" "$CPP_ROOT/cluster/results_${RUN_NAME}"
cd "$HERE"

echo ">>> RUN_NAME=$RUN_NAME  N=$n  NSEEDS=$NSEEDS  TIMEOUT=${TIMEOUT_OVERRIDE:-<csv>}"

sbid=$(
    export PROJECT_ROOT="$CPP_ROOT"
    export DATASETS_FILE DATA_BASE RUN_NAME NSEEDS TIMEOUT_OVERRIDE REGIME
    sbatch --account="$ACCOUNT" --array="1-${n}${ARRAY_LIMIT}" \
        --mem="$MEM" --cpus-per-task=1 --time="$TIME_LIMIT" \
        --output="logs/${RUN_NAME}_%A_%a.out" \
        --error="logs/${RUN_NAME}_%A_%a.err" \
        --export=ALL \
        job_boolmaxls_dataset.sh | awk '{print $NF}'
)
echo ">>> Job soumis : $sbid  (results_${RUN_NAME}/)"
