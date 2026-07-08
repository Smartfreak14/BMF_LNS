#!/bin/bash
# Compare opt4-wls contre baseline MaxSAT anytime sur une liste de datasets.
#
# Usage :
#   bench_maxsat_batch.sh <datasets.csv> <data_base> <timeout_s> <out_dir>
#
#   datasets.csv : filename,k,timeout_per_option_seconds (comme cluster/datasets*.csv)
#   data_base    : ex. Real/SimplUni
#   timeout_s    : timeout applique a CHAQUE solveur baseline
#   out_dir      : dossier de sortie ; chaque dataset produit :
#                    - <base>.wcnf
#                    - <base>_bench.csv (comparaison des solveurs sur cette instance)

set -euo pipefail

DATASETS="${1:?Usage: bench_maxsat_batch.sh <datasets.csv> <data_base> <timeout_s> <out_dir>}"
DATA_BASE="${2:?}"
TIMEOUT="${3:?}"
OUT_DIR="${4:?}"

HERE="$(cd "$(dirname "$0")" && pwd)"
CPP_ROOT="$(cd "$HERE/.." && pwd)"
RBAC_BIN="$CPP_ROOT/build/rbac_main"
DATA_DIR="$CPP_ROOT/data/$DATA_BASE"

[ -x "$RBAC_BIN" ] || { echo "ERREUR : $RBAC_BIN absent. Compiler d'abord."; exit 1; }
[ -d "$DATA_DIR" ] || { echo "ERREUR : $DATA_DIR absent."; exit 1; }

mkdir -p "$OUT_DIR"

# Lignes non commentees et non vides, en ignorant l'entete.
awk 'NR==1 {next} /^#/ {next} NF' "$DATASETS" | while IFS=',' read -r FILE K TOUT; do
    FILE="$(echo "$FILE" | tr -d '[:space:]')"
    K="$(echo "$K" | tr -d '[:space:]')"
    [ -f "$DATA_DIR/$FILE" ] || { echo "[skip] $FILE : absent"; continue; }

    BASE="${FILE%.csv}"
    WCNF="$OUT_DIR/${BASE}.wcnf"
    CSV="$OUT_DIR/${BASE}_bench.csv"

    echo ""
    echo "=================================================="
    echo " $FILE (k=$K, timeout=${TIMEOUT}s)"
    echo "=================================================="

    "$RBAC_BIN" --dump-wcnf "$DATA_DIR/$FILE" "$K" "$WCNF"
    bash "$HERE/bench_maxsat.sh" "$WCNF" "$TIMEOUT" "$CSV"
done

echo ""
echo ">>> Resultats individuels : $OUT_DIR/*_bench.csv"
