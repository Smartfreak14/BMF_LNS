#!/bin/bash
# Lance ASSO EN PARALLELE sur 31 datasets Real avec K = floor(K_papier / 2)
# (datasets_table2.csv). Protocole Table 2 du papier Florent.
#
# Meme que submit_asso_real.sh mais avec K/2 au lieu de K entier.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

DATASETS_FILE="${DATASETS_FILE:-$HERE/datasets_table2.csv}"
DATA_BASE="Real"
NSEEDS="${NSEEDS:-3}"
TIMEOUT_OVERRIDE="${TIMEOUT_OVERRIDE:-}"
REGIME="${REGIME:-raw}"
RUN_NAME="${RUN_NAME:-asso_table2}"

ACCOUNT="${ACCOUNT:-def-florenta}"
MEM="${MEM:-16G}"
TIME_LIMIT="${TIME_LIMIT:-12:00:00}"
ARRAY_LIMIT="${ARRAY_LIMIT:-}"

if [ -d "$HERE/../scripts/baselines" ] && [ -d "$HERE/../data" ]; then
    CPP_ROOT="$(cd "$HERE/.." && pwd)"
else
    echo "ERREUR : impossible de trouver scripts/ et data/ depuis $HERE" >&2
    exit 1
fi

[ -f "$CPP_ROOT/scripts/baselines/asso.py" ]        || { echo "ERREUR : asso.py absent"; exit 1; }
[ -f "$DATASETS_FILE" ]                              || { echo "ERREUR : $DATASETS_FILE absent"; exit 1; }
[ -d "$CPP_ROOT/data/$DATA_BASE" ]                   || { echo "ERREUR : data/$DATA_BASE absent"; exit 1; }
[ -f "$HERE/job_asso_dataset.sh" ]                   || { echo "ERREUR : job_asso_dataset.sh absent"; exit 1; }

n=$(awk 'NR==1 {next} /^#/ {next} NF {x++} END {print x+0}' "$DATASETS_FILE")
[ "$n" -gt 0 ] || { echo "ERREUR : $DATASETS_FILE vide"; exit 1; }

mkdir -p "$CPP_ROOT/cluster/logs" "$CPP_ROOT/cluster/results_${RUN_NAME}"
cd "$HERE"

array_spec="1-${n}${ARRAY_LIMIT}"

echo ">>> RUN_NAME         = $RUN_NAME"
echo ">>> DATASETS_FILE    = $DATASETS_FILE  (K/2)"
echo ">>> DATA_BASE        = $DATA_BASE"
echo ">>> N datasets       = $n"
echo ">>> Array spec       = $array_spec"
echo ">>> NSEEDS           = $NSEEDS"
echo ">>> TIMEOUT_OVERRIDE = ${TIMEOUT_OVERRIDE:-<from csv>}"
echo ">>> Time / task      = $TIME_LIMIT"
echo ">>> Mem / task       = $MEM"
echo ""

sbid=$(
    export PROJECT_ROOT="$CPP_ROOT"
    export DATASETS_FILE DATA_BASE RUN_NAME NSEEDS TIMEOUT_OVERRIDE REGIME
    sbatch \
        --account="$ACCOUNT" \
        --array="$array_spec" \
        --mem="$MEM" \
        --cpus-per-task=1 \
        --time="$TIME_LIMIT" \
        --output="logs/${RUN_NAME}_%A_%a.out" \
        --error="logs/${RUN_NAME}_%A_%a.err" \
        --export=ALL \
        job_asso_dataset.sh | awk '{print $NF}'
)

echo ">>> Job soumis    : $sbid"
echo ">>> Resultats     : $CPP_ROOT/cluster/results_${RUN_NAME}/<dataset>_asso_s<N>.csv"
echo ">>> Logs          : $CPP_ROOT/cluster/logs/${RUN_NAME}_${sbid}_<task>.{out,err}"
echo ""
echo ">>> Fin : python $CPP_ROOT/cluster/aggregate_table2.py \\"
echo "                 $CPP_ROOT/cluster/results_${RUN_NAME}/ -o ${RUN_NAME}_all.csv"
