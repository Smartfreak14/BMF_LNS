#!/bin/bash
# Ablation : 31 matrices brutes (data/Real/) a k = floor(k_real/2) x 5 configurations x NSEEDS graines.
# Configurations : greedy / no-weighting / wls-no-double-no-lns / no-lns / full.
# Timeout par config = colonne 3 de datasets_real_khalf.csv (1800s).
# Budget par tache : 4 configs x 1 seed x 1800s = 2h + verif.py (~15 min max sur advertisement).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

# EXPORTS AVANT source (sinon _submit_common.sh met TIME_LIMIT=00:30:00).
export JOB_SCRIPT="job_ablation_array.sh"
export NSEEDS="${NSEEDS:-1}"
export REGIME="${REGIME:-raw}"
export DATA_BASE="${DATA_BASE:-Real}"
export MEM_PER_TASK="${MEM_PER_TASK:-16G}"
export TIME_LIMIT="${TIME_LIMIT:-03:00:00}"

source "$HERE/_submit_common.sh"

submit_array \
    "$HERE/datasets_real_khalf.csv" \
    "$DATA_BASE" \
    "${RUN_NAME:-ablation_real_khalf}"
