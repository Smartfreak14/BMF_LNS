#!/bin/bash
# Baseline MaxSAT anytime sur les datasets racines data/Real/.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export JOB_SCRIPT="job_maxsat_array.sh"
export MAXSAT_TIMEOUT="${MAXSAT_TIMEOUT:-60}"

source "$HERE/_submit_common.sh"

submit_array \
    "$HERE/datasets_real.csv" \
    "Real" \
    "maxsat_real"
