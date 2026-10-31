#!/bin/bash
# 1 tache SLURM = 1 dataset x NSEEDS seeds pour Kolomvakis.
# Utilise par submit_kolomvakis_parallel.sh en mode array.
#
# Env attendus (via --export=ALL) :
#   PROJECT_ROOT      : racine du projet
#   DATASETS_FILE     : CSV filename,k,timeout
#   DATA_BASE         : "Real" (relatif a data/)
#   RUN_NAME          : suffixe pour results_<RUN_NAME>/
#   NSEEDS            : nombre de graines (defaut 3)
#   TIMEOUT_OVERRIDE  : timeout en secondes (defaut : utilise csv)
#   REGIME            : etiquette 'raw'/'simpluni'/etc

#SBATCH --account=def-florenta
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/kolomvakis_%A_%a.out
#SBATCH --error=logs/kolomvakis_%A_%a.err

set -euo pipefail

if [ -z "${PROJECT_ROOT:-}" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_ROOT="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
fi
: "${PROJECT_ROOT:?PROJECT_ROOT non defini}"

CPP_ROOT="$PROJECT_ROOT"
DATA_BASE="${DATA_BASE:-Real}"
RUN_NAME="${RUN_NAME:-kolomvakis}"
NSEEDS="${NSEEDS:-3}"
TIMEOUT_OVERRIDE="${TIMEOUT_OVERRIDE:-}"
REGIME="${REGIME:-raw}"

DATASETS_FILE="${DATASETS_FILE:-$CPP_ROOT/cluster/datasets_real.csv}"
DATA_DIR="$CPP_ROOT/data/$DATA_BASE"
RESULTS_DIR="$CPP_ROOT/cluster/results_${RUN_NAME}"
LOGS_DIR="$CPP_ROOT/cluster/logs"

mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

# Charger MATLAB + Gurobi (versions fixes pour eviter le prompt interactif)
module --force purge
module load StdEnv/2023 matlab/2026a.3 gurobi/11.0.0

# Recuperer la ligne du dataset a traiter
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
echo " EBROOTGUROBI = $EBROOTGUROBI"
echo "============================================================"

# Lancer MATLAB pour ce dataset, boucle interne sur NSEEDS
# ATTENTION : ici $EBROOTGUROBI est EVALUE cote shell (pas de \)
# car on est deja dans le job SLURM et le module gurobi est charge.
for SEED in $(seq 1 $NSEEDS); do
    OUT_CSV="$RESULTS_DIR/${BASE}_kolomvakis_s${SEED}.csv"
    TRACE="$RESULTS_DIR/${BASE}_kolomvakis_s${SEED}.trace"
    LOG="$RESULTS_DIR/${BASE}_kolomvakis_s${SEED}.log"
    echo ">>> seed=$SEED"
    # ATTENTION : pas de '...' MATLAB dans -batch, tout sur une ligne logique.
    # Le \ shell concatene deja tout en une commande.
    matlab -batch "addpath('$EBROOTGUROBI/matlab'); addpath('$CPP_ROOT/scripts/baselines'); run_kolomvakis('$CSV_PATH', $K, $TOUT, '$OUT_CSV', '$TRACE', $SEED, '$REGIME')" \
        > "$LOG" 2>&1 || echo "[warn] seed $SEED failed (voir $LOG)"
done

echo "Finished : $(date -Iseconds)"
