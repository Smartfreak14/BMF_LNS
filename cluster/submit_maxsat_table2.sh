#!/bin/bash
# Lance NuWLS (via bench_maxsat.sh) sur les 31 datasets Real avec K/2.
# Meme structure que submit_maxsat_real.sh mais utilise datasets_table2.csv.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

export JOB_SCRIPT="job_maxsat_array.sh"
export MAXSAT_TIMEOUT="${MAXSAT_TIMEOUT:-60}"
export DATASETS_FILE="${DATASETS_FILE:-$HERE/datasets_table2.csv}"
export RUN_NAME="${RUN_NAME:-maxsat_table2}"

source "$HERE/_submit_common.sh"

submit_array \
    "$DATASETS_FILE" \
    "Real" \
    "$RUN_NAME"
