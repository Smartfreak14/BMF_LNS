#!/bin/bash
# Lance BoolMaxLS (BoolMaxLS = notre methode) en PARALLELE sur 31 datasets Real
# avec K entier (datasets_real.csv).
# 1 tache SLURM = 1 dataset x NSEEDS graines.
#
# Env :
#   NSEEDS           : graines (defaut 3)
#   TIMEOUT_OVERRIDE : timeout secondes (defaut : csv)
#   REGIME           : etiquette regime (defaut 'raw')
#   RUN_NAME         : suffixe results_<RUN_NAME>/ (defaut boolmaxls_real)
#   SEED_PER_TASK    : 1 = une tache par (dataset, graine) ; INIT_ALGO, BOOLMAXLS_ARGS : voir job_boolmaxls_dataset.sh

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

DATASETS_FILE="${DATASETS_FILE:-$HERE/datasets_real.csv}"
DATA_BASE="${DATA_BASE:-Real}"
NSEEDS="${NSEEDS:-3}"
TIMEOUT_OVERRIDE="${TIMEOUT_OVERRIDE:-}"
REGIME="${REGIME:-raw}"
RUN_NAME="${RUN_NAME:-boolmaxls_real}"
SEED_PER_TASK="${SEED_PER_TASK:-0}"
INIT_ALGO="${INIT_ALGO:-}"
BOOLMAXLS_ARGS="${BOOLMAXLS_ARGS:-}"

ACCOUNT="${ACCOUNT:-def-florenta}"
MEM="${MEM:-16G}"
TIME_LIMIT="${TIME_LIMIT:-12:00:00}"
ARRAY_LIMIT="${ARRAY_LIMIT:-}"

if [ -d "$HERE/../src" ] && [ -d "$HERE/../data" ]; then
    CPP_ROOT="$(cd "$HERE/.." && pwd)"
else
    echo "ERREUR : structure introuvable depuis $HERE" >&2
    exit 1
fi

[ -x "$CPP_ROOT/build/BoolMaxLS" ] || { echo "ERREUR : build/BoolMaxLS absent"; exit 1; }
[ -f "$DATASETS_FILE" ]             || { echo "ERREUR : $DATASETS_FILE absent"; exit 1; }
[ -d "$CPP_ROOT/data/$DATA_BASE" ]   || { echo "ERREUR : data/$DATA_BASE absent"; exit 1; }
[ -f "$HERE/job_boolmaxls_dataset.sh" ] || { echo "ERREUR : job_boolmaxls_dataset.sh absent"; exit 1; }

n=$(awk 'NR==1 {next} /^#/ {next} NF {x++} END {print x+0}' "$DATASETS_FILE")
[ "$SEED_PER_TASK" = 1 ] && n=$(( n * NSEEDS ))
[ -z "$INIT_ALGO" ] || [ -x "$CPP_ROOT/build/FromBelowBaselines" ] || { echo "ERREUR : build/FromBelowBaselines absent"; exit 1; }

mkdir -p "$CPP_ROOT/cluster/logs" "$CPP_ROOT/cluster/results_${RUN_NAME}"
cd "$HERE"

array_spec="1-${n}${ARRAY_LIMIT}"

echo ">>> RUN_NAME=$RUN_NAME  N=$n  NSEEDS=$NSEEDS  TIMEOUT=${TIMEOUT_OVERRIDE:-<csv>}  Time=$TIME_LIMIT"

sbid=$(
    export PROJECT_ROOT="$CPP_ROOT"
    export DATASETS_FILE DATA_BASE RUN_NAME NSEEDS TIMEOUT_OVERRIDE REGIME SEED_PER_TASK INIT_ALGO BOOLMAXLS_ARGS
    sbatch \
        --account="$ACCOUNT" --array="$array_spec" \
        --mem="$MEM" --cpus-per-task=1 --time="$TIME_LIMIT" \
        --output="logs/${RUN_NAME}_%A_%a.out" \
        --error="logs/${RUN_NAME}_%A_%a.err" \
        --export=ALL \
        job_boolmaxls_dataset.sh | awk '{print $NF}'
)

echo ">>> Job soumis : $sbid  (results_${RUN_NAME}/)"
echo "$RUN_NAME $sbid $(date '+%Y-%m-%d %H:%M')" >> "$CPP_ROOT/cluster/submitted_jobs.log"
