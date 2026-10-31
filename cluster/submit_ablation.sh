#!/bin/bash
# Ablation study : 31 datasets Real/SimplUni x 5 configurations x NSEEDS graines.
# Configurations : greedy / no-weighting / wls-no-double-no-lns / no-lns / full.
# Timeout par config = colonne 3 du CSV (1800s = 30 min par defaut).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

# EXPORTS AVANT source (sinon _submit_common.sh met TIME_LIMIT=00:30:00).
export JOB_SCRIPT="job_ablation_array.sh"
export NSEEDS="${NSEEDS:-5}"                # avec 5 configs x 1800s -> 27000s = 7.5h
export REGIME="${REGIME:-simpluni}"
export DATA_BASE="${DATA_BASE:-Real/SimplUni}"
export MEM_PER_TASK="${MEM_PER_TASK:-16G}"
export TIME_LIMIT="${TIME_LIMIT:-00:40:00}"

source "$HERE/_submit_common.sh"

submit_array \
    "$HERE/datasets_ablation.csv" \
    "$DATA_BASE" \
    "ablation_simpluni"
