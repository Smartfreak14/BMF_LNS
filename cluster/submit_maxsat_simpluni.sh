#!/bin/bash
# Baseline MaxSAT anytime sur tous les datasets SimplUni.
# Timeout par solveur : MAXSAT_TIMEOUT (defaut 60s).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export JOB_SCRIPT="job_maxsat_array.sh"
export MAXSAT_TIMEOUT="${MAXSAT_TIMEOUT:-60}"

source "$HERE/_submit_common.sh"

submit_array \
    "$HERE/datasets.csv" \
    "Real/SimplUni" \
    "maxsat_simpluni"
