#!/bin/bash
# Lance Kolomvakis (BoolMF+CombineBoolMFs) en batch sur les 31 datasets Real.
# Une seule tache SLURM qui parcourt les datasets en sequence.
#
# Prerequis :
#   - Gurobi + MATLAB fonctionnels : gurobi_cl --license (pas d'erreur)
#   - matlab -batch "gurobi_setup" execute au moins une fois
#   - Sources Kolomvakis dans third_party/BooleanMF/
#   - Wrappers scripts/baselines/run_kolomvakis.m et run_kolomvakis_batch.m
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

DATASETS_FILE="${DATASETS_FILE:-$HERE/datasets_real.csv}"
DATA_DIR_REL="${DATA_DIR_REL:-Real}"
NSEEDS="${NSEEDS:-5}"
TIMEOUT_OVERRIDE="${TIMEOUT_OVERRIDE:-}"       # secondes, vide = pas d'override
REGIME="${REGIME:-raw}"
RUN_NAME="${RUN_NAME:-kolomvakis_batch}"

ACCOUNT="${ACCOUNT:-def-florenta}"
MEM="${MEM:-32G}"
TIME_LIMIT="${TIME_LIMIT:-23:00:00}"

# Racine du projet (contient scripts/, data/, third_party/)
if [ -d "$HERE/../scripts/baselines" ] && [ -d "$HERE/../data" ]; then
    CPP_ROOT="$(cd "$HERE/.." && pwd)"
else
    echo "ERREUR : impossible de trouver scripts/baselines/ et data/ depuis $HERE" >&2
    exit 1
fi

# Sanity checks
[ -f "$CPP_ROOT/scripts/baselines/run_kolomvakis.m" ]        || { echo "ERREUR : run_kolomvakis.m absent"; exit 1; }
[ -f "$CPP_ROOT/scripts/baselines/run_kolomvakis_batch.m" ]  || { echo "ERREUR : run_kolomvakis_batch.m absent"; exit 1; }
[ -d "$CPP_ROOT/third_party/BooleanMF" ]                     || { echo "ERREUR : third_party/BooleanMF/ absent"; exit 1; }
[ -f "$DATASETS_FILE" ]                                       || { echo "ERREUR : $DATASETS_FILE absent"; exit 1; }
[ -d "$CPP_ROOT/data/$DATA_DIR_REL" ]                         || { echo "ERREUR : data/$DATA_DIR_REL absent"; exit 1; }

OUTPUT_DIR="$CPP_ROOT/cluster/results_${RUN_NAME}"
mkdir -p "$OUTPUT_DIR" "$CPP_ROOT/cluster/logs"

# Argument MATLAB timeout : [] si vide, sinon le nombre
if [ -z "$TIMEOUT_OVERRIDE" ]; then
    TO_ARG="[]"
else
    TO_ARG="$TIMEOUT_OVERRIDE"
fi

echo ">>> Kolomvakis batch : $DATASETS_FILE"
echo "    data_dir   = $CPP_ROOT/data/$DATA_DIR_REL"
echo "    output_dir = $OUTPUT_DIR"
echo "    nseeds     = $NSEEDS"
echo "    timeout    = ${TIMEOUT_OVERRIDE:-<from csv>}"
echo "    regime     = $REGIME"
echo ""

# Note : on charge matlab + gurobi dans le wrap.
# Le module gurobi definit $EBROOTGUROBI vers l'installation.
sbid=$(sbatch \
    --account="$ACCOUNT" \
    --mem="$MEM" \
    --cpus-per-task=1 \
    --time="$TIME_LIMIT" \
    --output="$CPP_ROOT/cluster/logs/${RUN_NAME}_%j.out" \
    --error="$CPP_ROOT/cluster/logs/${RUN_NAME}_%j.err" \
    --wrap="\
module load StdEnv/2023 matlab/2026a.3 gurobi/11.0.0 && \
matlab -batch \"\
    addpath('\$EBROOTGUROBI/matlab'); \
    addpath('$CPP_ROOT/scripts/baselines'); \
    run_kolomvakis_batch('$DATASETS_FILE', '$CPP_ROOT/data/$DATA_DIR_REL', '$OUTPUT_DIR', $NSEEDS, $TO_ARG, '$REGIME')\
\"" \
    | awk '{print $NF}')

echo ">>> Job soumis : $sbid"
echo ">>> Log        : $CPP_ROOT/cluster/logs/${RUN_NAME}_${sbid}.out"
echo ">>> Resultats  : $OUTPUT_DIR/<dataset>_kolomvakis_s<N>.csv"
echo ""
echo ">>> Suivi : tail -f $CPP_ROOT/cluster/logs/${RUN_NAME}_${sbid}.out"
echo ">>> Une fois fini : python $CPP_ROOT/cluster/aggregate_table2.py $OUTPUT_DIR -o kolomvakis_all.csv"
