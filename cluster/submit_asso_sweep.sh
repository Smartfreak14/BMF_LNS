#!/bin/bash
# ASSO reviewer-proof : lance ASSO en balayant explicitement tau sur
# {0.3, 0.5, 0.7, 0.9} (plage calibree par Kolomvakis, etendue a 0.9
# pour les matrices denses). Budget partage equitablement : chaque tau
# recoit timeout/4 secondes -> ASSO consomme le meme budget total que
# les autres baselines.
#
# Pour chaque (dataset, seed), asso.py produit :
#   - <base>_asso_s<N>.csv    : best-of tau (schema unifie)
#   - <base>_asso_s<N>.trace  : trace time_ms err (best current)
#   - <base>_asso_s<N>.tau.csv: erreur par tau (audit reviewer)
#
# Env :
#   TAU_LIST         : liste virgule-separee (defaut 0.3,0.5,0.7,0.9)
#   NSEEDS           : graines / dataset (defaut 3)
#   TIMEOUT_OVERRIDE : timeout total en s (defaut 1800s, meme budget que BoolMaxLS)
#   DATASETS_FILE    : CSV filename,k,timeout (defaut datasets_real.csv)
#   DATA_BASE        : sous-dossier data/ (defaut Real)
#   REGIME           : etiquette (defaut raw)

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

export TAU_LIST="${TAU_LIST:-0.3,0.5,0.7,0.9}"
export DATASETS_FILE="${DATASETS_FILE:-$HERE/datasets_real.csv}"
export DATA_BASE="${DATA_BASE:-Real}"
export NSEEDS="${NSEEDS:-3}"
export TIMEOUT_OVERRIDE="${TIMEOUT_OVERRIDE:-1800}"
export REGIME="${REGIME:-raw}"
export RUN_NAME="${RUN_NAME:-asso_sweep}"

export ACCOUNT="${ACCOUNT:-def-florenta}"
export MEM="${MEM:-16G}"
export TIME_LIMIT="${TIME_LIMIT:-12:00:00}"
export ARRAY_LIMIT="${ARRAY_LIMIT:-}"

CPP_ROOT="$(cd "$HERE/.." && pwd)"

[ -f "$CPP_ROOT/scripts/baselines/asso.py" ] || { echo "ERREUR : asso.py absent"; exit 1; }
[ -f "$DATASETS_FILE" ]                       || { echo "ERREUR : $DATASETS_FILE absent"; exit 1; }
[ -d "$CPP_ROOT/data/$DATA_BASE" ]            || { echo "ERREUR : data/$DATA_BASE absent"; exit 1; }
[ -f "$HERE/job_asso_dataset.sh" ]            || { echo "ERREUR : job_asso_dataset.sh absent"; exit 1; }

n=$(awk 'NR==1 {next} /^#/ {next} NF {x++} END {print x+0}' "$DATASETS_FILE")
[ "$n" -gt 0 ] || { echo "ERREUR : $DATASETS_FILE vide"; exit 1; }

n_tau=$(awk -F, '{print NF}' <<< "$TAU_LIST")
per_tau=$(awk "BEGIN {printf \"%.0f\", $TIMEOUT_OVERRIDE / $n_tau}")

mkdir -p "$CPP_ROOT/cluster/logs" "$CPP_ROOT/cluster/results_${RUN_NAME}"
cd "$HERE"

array_spec="1-${n}${ARRAY_LIMIT}"

echo ">>> RUN_NAME          = $RUN_NAME"
echo ">>> TAU_LIST          = $TAU_LIST  (n_tau=$n_tau)"
echo ">>> TIMEOUT total     = ${TIMEOUT_OVERRIDE}s"
echo ">>> Budget par tau    = ${per_tau}s"
echo ">>> DATASETS_FILE     = $DATASETS_FILE ($n datasets)"
echo ">>> NSEEDS            = $NSEEDS"
echo ">>> Time / task SLURM = $TIME_LIMIT"
echo ""

sbid=$(
    export PROJECT_ROOT="$CPP_ROOT"
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
echo ">>> Resultats     : $CPP_ROOT/cluster/results_${RUN_NAME}/"
echo ">>> Audit tau     : results_${RUN_NAME}/<dataset>_asso_s<N>.tau.csv"
echo ""
echo ">>> Suivi : squeue -u \$USER | grep ${RUN_NAME:0:8}"
echo ">>> Fin   : python $CPP_ROOT/cluster/aggregate_table2.py \\"
echo "                  $CPP_ROOT/cluster/results_${RUN_NAME}/ -o ${RUN_NAME}_all.csv"
