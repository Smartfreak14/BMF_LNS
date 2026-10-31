#!/bin/bash
# Sourcé par tous les submit_*.sh. Définit submit_array() qui pousse vers
# SLURM avec les bons env vars + ressources.
#
# Variables attendues côté appelant :
#   PROJECT_ROOT    (auto si non set)
#   ACCOUNT         (defaut def-florenta)
#   MEM_PER_TASK    (defaut 8G)
#   CPUS_PER_TASK   (defaut 1)
#   TIME_LIMIT      (defaut 23:00:00)
#   ARRAY_LIMIT     (ex. "%8" pour 8 max en parallele)
#   SEED_PER_TASK   (1 = N x NSEEDS taches, une par dataset et graine)
#
# Usage : submit_array <datasets_file> <data_base> <run_name>
#   data_base est relatif a cpp/data/ (ex. "Real/SimplUni", "Synthetic").

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
if [ -f "$PROJECT_ROOT/CMakeLists.txt" ]; then
    CPP_ROOT="$PROJECT_ROOT"
elif [ -f "$PROJECT_ROOT/BoolMaxLStenance/cpp/CMakeLists.txt" ]; then
    CPP_ROOT="$PROJECT_ROOT/BoolMaxLStenance/cpp"
else
    echo "ERREUR : pas de CMakeLists.txt sous $PROJECT_ROOT"
    exit 1
fi
CLUSTER_DIR="$CPP_ROOT/cluster"

ACCOUNT="${ACCOUNT:-def-florenta}"
MEM_PER_TASK="${MEM_PER_TASK:-32G}"
CPUS_PER_TASK="${CPUS_PER_TASK:-1}"
TIME_LIMIT="${TIME_LIMIT:-00:30:00}"
ARRAY_LIMIT="${ARRAY_LIMIT:-}"

submit_array() {
    local datasets_file="$1"
    local data_base="$2"
    local run_name="$3"

    [ -x "$CPP_ROOT/build/BoolMaxLS" ] || { echo "ERREUR : $CPP_ROOT/build/BoolMaxLS absent. bash setup_env.sh d'abord."; exit 1; }
    [ -f "$datasets_file" ]            || { echo "ERREUR : $datasets_file introuvable.";                                    exit 1; }
    [ -d "$CPP_ROOT/data/$data_base" ] || { echo "ERREUR : dossier de donnees $CPP_ROOT/data/$data_base introuvable.";      exit 1; }

    local n
    n=$(awk 'NR==1 {next} /^#/ {next} NF {x++} END {print x+0}' "$datasets_file")
    if [ "$n" -le 0 ]; then
        echo "ERREUR : $datasets_file vide (N=$n)."
        exit 1
    fi
    # SEED_PER_TASK=1 : une tache par (dataset, graine).
    [ "${SEED_PER_TASK:-0}" != 1 ] || n=$(( n * NSEEDS ))

    mkdir -p "$CLUSTER_DIR/logs" "$CLUSTER_DIR/results_${run_name}"
    cd "$CLUSTER_DIR"

    local array_spec="1-${n}${ARRAY_LIMIT}"

    echo ">>> RUN_NAME      = $run_name"
    echo ">>> PROJECT_ROOT  = $PROJECT_ROOT"
    echo ">>> DATASETS_FILE = $datasets_file"
    echo ">>> DATA_BASE     = $data_base"
    echo ">>> N datasets    = $n"
    echo ">>> Array spec    = $array_spec"
    echo ">>> Account       = $ACCOUNT"
    echo ">>> Mem / task    = $MEM_PER_TASK"
    echo ">>> CPUs / task   = $CPUS_PER_TASK"
    echo ">>> Time / task   = $TIME_LIMIT"
    echo ""

    local job_script="${JOB_SCRIPT:-job_array.sh}"

    local sbid
    sbid=$(
        export PROJECT_ROOT="$PROJECT_ROOT"
        export DATASETS_FILE="$datasets_file"
        export DATA_BASE="$data_base"
        export RUN_NAME="$run_name"
        export MAXSAT_TIMEOUT="${MAXSAT_TIMEOUT:-}"
        sbatch \
            --account="$ACCOUNT" \
            --array="$array_spec" \
            --mem="$MEM_PER_TASK" \
            --cpus-per-task="$CPUS_PER_TASK" \
            --time="$TIME_LIMIT" \
            --output="logs/${run_name}_%A_%a.out" \
            --error="logs/${run_name}_%A_%a.err" \
            --export=ALL \
            "$job_script" | awk '{print $NF}'
    )

    echo ">>> Job soumis  : $sbid"
    echo "$run_name $sbid $(date '+%Y-%m-%d %H:%M')" >> "$CLUSTER_DIR/submitted_jobs.log"
    echo ">>> Resultats   : $CLUSTER_DIR/results_${run_name}/<dataset>.csv"
    echo ">>> Logs        : $CLUSTER_DIR/logs/${run_name}_${sbid}_<task>.{out,err}"
    echo ""
    echo ">>> Apres la fin : python aggregate.py results_${run_name}/ -o results_${run_name}_all.csv"
}
