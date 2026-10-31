#!/bin/bash
# 1 tache SLURM = 1 dataset pour une baseline from-below (GreConD, GreEss, MDLGreConD).
# Nombre de facteurs = k de la liste (raw k ou raw floor(k/2)), ou, avec EXACT=1,
# sans limite : l'algorithme va jusqu'a la factorisation exacte (Simpli∀, Simpli∃₀).
# Chaque sortie est verifiee par src/verif.py.
#
# Env :
#   PROJECT_ROOT, DATASETS_FILE (filename,k,timeout), DATA_BASE, RUN_NAME
#   ALGO          : grecond | greess | mdlgrecond
#   EXACT         : 1 = factorisation exacte (defaut 0)
#   BUDGET        : limite de temps de l'algorithme en secondes (defaut 1800)
#   MISSING_MODES : traitements des cellules inconnues (defaut "free zero one") ;
#                   sans cellule inconnue les trois sont identiques, seul free est lance

#SBATCH --account=def-florenta
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

if [ -z "${PROJECT_ROOT:-}" ] && [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    PROJECT_ROOT="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
fi
: "${PROJECT_ROOT:?PROJECT_ROOT non defini}"
: "${ALGO:?ALGO non defini (grecond | greess | mdlgrecond)}"

CPP_ROOT="$PROJECT_ROOT"
DATA_BASE="${DATA_BASE:-Real}"
RUN_NAME="${RUN_NAME:-${ALGO}_raw_k}"
EXACT="${EXACT:-0}"
BUDGET="${BUDGET:-1800}"
MISSING_MODES="${MISSING_MODES:-free zero one}"

DATASETS_FILE="${DATASETS_FILE:-$CPP_ROOT/cluster/datasets_real.csv}"
BIN="$CPP_ROOT/build/FromBelowBaselines"
DATA_DIR="$CPP_ROOT/data/$DATA_BASE"
RESULTS_DIR="$CPP_ROOT/cluster/results_${RUN_NAME}"
mkdir -p "$RESULTS_DIR" "$CPP_ROOT/cluster/logs"

module --force purge
module load StdEnv/2023 gcc/12.3 python/3.11

TASK_ID="${SLURM_ARRAY_TASK_ID:-}"
[ -n "$TASK_ID" ] || { echo "ERREUR : SLURM_ARRAY_TASK_ID absent"; exit 1; }
LINE=$(grep -v '^#' "$DATASETS_FILE" | tail -n +2 | sed -n "${TASK_ID}p")
[ -n "$LINE" ] || { echo "ERREUR : ligne $TASK_ID absente de $DATASETS_FILE"; exit 1; }

IFS=',' read -r FILE K _IGNORED <<< "$LINE"
FILE="$(echo "$FILE" | tr -d '[:space:]')"
K="$(echo "$K" | tr -d '[:space:]')"
CSV_PATH="$DATA_DIR/$FILE"
[ -f "$CSV_PATH" ] || { echo "ERREUR : $CSV_PATH absent"; exit 1; }
BASE="${FILE%.csv}"
CUT="$K"
[ "$EXACT" = 1 ] && CUT=1000000

# Cellule inconnue = champ vide ou "?" (meme convention que CSVMatrixLoader).
HAS_UNKNOWN=$(awk -F',' '{for (i = 1; i <= NF; i++) { g = $i; gsub(/[ \t\r]/, "", g); if (g == "" || g == "?") { print 1; exit } } }' "$CSV_PATH")

echo "=== $ALGO  $DATA_BASE/$FILE  facteurs<=$CUT  budget=${BUDGET}s  inconnues=${HAS_UNKNOWN:-0} ==="

for MISS in $MISSING_MODES; do
    if [ "$MISS" != "free" ] && [ -z "$HAS_UNKNOWN" ]; then continue; fi
    STEM="$RESULTS_DIR/${BASE}_${ALGO}_${MISS}"
    "$BIN" "$ALGO" "$CSV_PATH" "$CUT" --missing "$MISS" --time-limit "$BUDGET" \
        --out-dir "$RESULTS_DIR" > "$STEM.log" 2>&1 || echo "[warn] $ALGO $MISS : code $? (voir $STEM.log)"

    echo "k,verif_k,covered_zeros,uncovered_ones,error" > "${STEM}_verif.csv"
    A="${STEM}_k${CUT}_A.csv"; B="${STEM}_k${CUT}_B.csv"
    [ -f "$A" ] || continue
    if [ ! -s "$B" ]; then
        echo "$CUT,0,0,NA,NA" >> "${STEM}_verif.csv"   # aucun facteur (budget epuise avant le premier)
        continue
    fi
    OUT=$(cd "$CPP_ROOT/src" && python3 verif.py "$A" "$B" "$RESULTS_DIR/${BASE}_M.csv" 2>&1) || true
    echo "$OUT" > "${STEM}_verif.log"
    get() { echo "$OUT" | sed -n "s/^$1 = //p"; }
    echo "$CUT,$(get 'k'),$(get 'Covered 0-entries (hard violations)'),$(get 'Uncovered 1-entries (soft violations, err_M)'),$(get 'Reconstruction error')" >> "${STEM}_verif.csv"
    gzip -f "$A" "$B"
done
rm -f "$RESULTS_DIR/${BASE}_M.csv"

echo "Finished : $(date -Iseconds)"
