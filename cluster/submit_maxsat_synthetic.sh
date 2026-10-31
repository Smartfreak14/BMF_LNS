#!/bin/bash
# Baseline MaxSAT anytime sur tous les datasets Synthetic.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

export JOB_SCRIPT="job_maxsat_array.sh"
export MAXSAT_TIMEOUT="${MAXSAT_TIMEOUT:-60}"

source "$HERE/_submit_common.sh"

DATASETS_FILE="${DATASETS_FILE:-$HERE/datasets_synthetic.csv}"
SYNTHETIC_DIR="$CPP_ROOT/data/Synthetic"
TIMEOUT_PER_OPT="${TIMEOUT_PER_OPT:-60}"

bash "$HERE/build_synthetic_list.sh" "$SYNTHETIC_DIR" "$TIMEOUT_PER_OPT" > "$DATASETS_FILE"

submit_array \
    "$DATASETS_FILE" \
    "Synthetic" \
    "maxsat_synthetic"
