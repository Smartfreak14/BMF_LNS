#!/bin/bash
# Lance Kolomvakis EN PARALLELE : 1 tache SLURM par dataset via un array job.
# Contrairement a submit_kolomvakis_batch.sh (1 tache pour les 31 en sequence),
# celui-ci lance 31 taches independantes qui tournent sur des noeuds differents.
#
# Chaque tache : 1 dataset x NSEEDS seeds, timeout par run configurable.
#
# Env :
#   NSEEDS            : graines par dataset (defaut 3)
#   TIMEOUT_OVERRIDE  : timeout secondes par run (defaut : csv)
#   REGIME            : etiquette regime (defaut 'raw')
#   RUN_NAME          : suffixe results_<RUN_NAME>/ (defaut kolomvakis_parallel)
#   DATASETS_FILE     : CSV filename,k,timeout (defaut datasets_real.csv)
#   DATA_BASE         : sous-dossier de data/ (defaut Real)
#   ARRAY_LIMIT       : ex "%8" pour 8 taches en parallele max
#   TIME_LIMIT        : temps SLURM par tache (defaut 05:00:00 = 5h)
#   MEM               : memoire par tache (defaut 32G)

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

DATASETS_FILE="${DATASETS_FILE:-$HERE/datasets_real.csv}"
DATA_BASE="${DATA_BASE:-Real}"
NSEEDS="${NSEEDS:-3}"
TIMEOUT_OVERRIDE="${TIMEOUT_OVERRIDE:-}"
REGIME="${REGIME:-raw}"
RUN_NAME="${RUN_NAME:-kolomvakis_parallel}"

ACCOUNT="${ACCOUNT:-def-florenta}"
MEM="${MEM:-32G}"
TIME_LIMIT="${TIME_LIMIT:-05:00:00}"
ARRAY_LIMIT="${ARRAY_LIMIT:-}"

# Racine du projet (contient scripts/, data/, third_party/)
if [ -d "$HERE/../scripts/baselines" ] && [ -d "$HERE/../data" ]; then
    CPP_ROOT="$(cd "$HERE/.." && pwd)"
else
    echo "ERREUR : impossible de trouver scripts/ et data/ depuis $HERE" >&2
    exit 1
fi

# Sanity checks
[ -f "$CPP_ROOT/scripts/baselines/run_kolomvakis.m" ] || { echo "ERREUR : run_kolomvakis.m absent"; exit 1; }
[ -d "$CPP_ROOT/third_party/BooleanMF" ]              || { echo "ERREUR : third_party/BooleanMF/ absent"; exit 1; }
[ -f "$CPP_ROOT/third_party/BooleanMF/utils/vec.m" ]  || { echo "ERREUR : third_party/BooleanMF/utils/vec.m absent"; exit 1; }
[ -f "$DATASETS_FILE" ]                                || { echo "ERREUR : $DATASETS_FILE absent"; exit 1; }
[ -d "$CPP_ROOT/data/$DATA_BASE" ]                     || { echo "ERREUR : data/$DATA_BASE absent"; exit 1; }
[ -f "$HERE/job_kolomvakis_dataset.sh" ]               || { echo "ERREUR : job_kolomvakis_dataset.sh absent"; exit 1; }

# Nombre de datasets (lignes non-commentees, hors entete)
n=$(awk 'NR==1 {next} /^#/ {next} NF {x++} END {print x+0}' "$DATASETS_FILE")
[ "$n" -gt 0 ] || { echo "ERREUR : $DATASETS_FILE vide"; exit 1; }

mkdir -p "$CPP_ROOT/cluster/logs" "$CPP_ROOT/cluster/results_${RUN_NAME}"
cd "$HERE"

array_spec="1-${n}${ARRAY_LIMIT}"

echo ">>> RUN_NAME         = $RUN_NAME"
echo ">>> PROJECT_ROOT     = $CPP_ROOT"
echo ">>> DATASETS_FILE    = $DATASETS_FILE"
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
    export DATASETS_FILE
    export DATA_BASE
    export RUN_NAME
    export NSEEDS
    export TIMEOUT_OVERRIDE
    export REGIME
    sbatch \
        --account="$ACCOUNT" \
        --array="$array_spec" \
        --mem="$MEM" \
        --cpus-per-task=1 \
        --time="$TIME_LIMIT" \
        --output="logs/${RUN_NAME}_%A_%a.out" \
        --error="logs/${RUN_NAME}_%A_%a.err" \
        --export=ALL \
        job_kolomvakis_dataset.sh | awk '{print $NF}'
)

echo ">>> Job soumis    : $sbid"
echo ">>> Resultats     : $CPP_ROOT/cluster/results_${RUN_NAME}/<dataset>_kolomvakis_s<N>.csv"
echo ">>> Logs SLURM    : $CPP_ROOT/cluster/logs/${RUN_NAME}_${sbid}_<task>.{out,err}"
echo ">>> Logs MATLAB   : $CPP_ROOT/cluster/results_${RUN_NAME}/<dataset>_kolomvakis_s<N>.log"
echo ""
echo ">>> Suivi : squeue -u \$USER | grep ${RUN_NAME:0:8}"
echo ">>> Fin   : python $CPP_ROOT/cluster/aggregate_table2.py \\"
echo "                  $CPP_ROOT/cluster/results_${RUN_NAME}/ -o ${RUN_NAME}_all.csv"
