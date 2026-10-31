#!/bin/bash
# MASTER SCRIPT : lance les 24 campagnes SLURM de la Table 2.
#
#   4 methodes : boolmaxls, nuwls, kolomvakis, asso
#   3 budgets  : 30s, 300s (5 min), 1800s (30 min)
#   2 rangs    : K entier (datasets_real.csv), K/2 (datasets_table2.csv)
#
# TOTAL : 24 campagnes independantes. Chacune est un array job de 31 taches
# SLURM en parallele = 24 * 31 = 744 taches SLURM au total.
#
# Chaque tache execute NSEEDS graines. Convention de nommage des sorties :
#   results_<method>_<krank>_<budget>/
# ex : results_boolmaxls_real_30s/, results_nuwls_k2_1800s/, etc.
#
# Options :
#   --dry-run     : affiche les commandes sans les executer
#   --method M    : lance uniquement M (boolmaxls|nuwls|kolomvakis|asso)
#   --budget B    : lance uniquement B (30|300|1800)
#   --rank R      : lance uniquement R (real|k2)
#   --nseeds N    : graines (defaut 3)
#
# Exemples :
#   bash launch_table2_all.sh                          # tout
#   bash launch_table2_all.sh --dry-run                # apercoit
#   bash launch_table2_all.sh --method boolmaxls         # BoolMaxLS uniquement (6 campagnes)
#   bash launch_table2_all.sh --budget 300             # tous les 5 min (8 campagnes)
#   bash launch_table2_all.sh --method asso --rank k2  # ASSO K/2 (3 campagnes)

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

DRY_RUN=0
METHOD_FILTER=""
BUDGET_FILTER=""
RANK_FILTER=""
NSEEDS="${NSEEDS:-3}"

while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run)     DRY_RUN=1 ;;
        --method)      METHOD_FILTER="$2"; shift ;;
        --budget)      BUDGET_FILTER="$2"; shift ;;
        --rank)        RANK_FILTER="$2"; shift ;;
        --nseeds)      NSEEDS="$2"; shift ;;
        *)             echo "arg inconnu : $1"; exit 1 ;;
    esac
    shift
done

BUDGETS="30 300 1800"
RANKS="real k2"
METHODS="boolmaxls nuwls kolomvakis asso"

# SLURM --time par budget (large marge pour NSEEDS seeds sequentiels)
declare -A SLURM_TIME
SLURM_TIME[30]="00:30:00"       # 30s * 3 seeds = 90s + overhead
SLURM_TIME[300]="03:00:00"      # 300s * 3 = 900s + overhead
SLURM_TIME[1800]="12:00:00"     # 1800s * 3 = 5400s + overhead

launch() {
    local method="$1" rank="$2" budget="$3"
    local time_lim="${SLURM_TIME[$budget]}"
    local run_name="${method}_${rank}_${budget}s"

    # Fichier datasets selon le rank
    local dataset_file
    if [ "$rank" = "real" ]; then
        dataset_file="$HERE/datasets_real.csv"
    else
        dataset_file="$HERE/datasets_table2.csv"
    fi

    # Commande selon la methode
    local cmd
    case "$method" in
        boolmaxls)
            if [ "$rank" = "real" ]; then
                cmd="NSEEDS=$NSEEDS TIMEOUT_OVERRIDE=$budget REGIME=raw"
                cmd="$cmd RUN_NAME=$run_name TIME_LIMIT=$time_lim"
                cmd="$cmd bash submit_boolmaxls_real.sh"
            else
                cmd="NSEEDS=$NSEEDS TIMEOUT_OVERRIDE=$budget REGIME=raw"
                cmd="$cmd RUN_NAME=$run_name TIME_LIMIT=$time_lim"
                cmd="$cmd bash submit_boolmaxls_table2.sh"
            fi
            ;;
        nuwls)
            if [ "$rank" = "real" ]; then
                cmd="MAXSAT_TIMEOUT=$budget DATASETS_FILE=$dataset_file"
                cmd="$cmd RUN_NAME=$run_name TIME_LIMIT=$time_lim"
                cmd="$cmd bash submit_maxsat_real.sh"
            else
                cmd="MAXSAT_TIMEOUT=$budget DATASETS_FILE=$dataset_file"
                cmd="$cmd RUN_NAME=$run_name TIME_LIMIT=$time_lim"
                cmd="$cmd bash submit_maxsat_table2.sh"
            fi
            ;;
        kolomvakis)
            cmd="NSEEDS=$NSEEDS TIMEOUT_OVERRIDE=$budget REGIME=raw"
            cmd="$cmd DATASETS_FILE=$dataset_file"
            cmd="$cmd RUN_NAME=$run_name TIME_LIMIT=$time_lim"
            cmd="$cmd bash submit_kolomvakis_parallel.sh"
            ;;
        asso)
            if [ "$rank" = "real" ]; then
                cmd="NSEEDS=$NSEEDS TIMEOUT_OVERRIDE=$budget REGIME=raw"
                cmd="$cmd RUN_NAME=$run_name TIME_LIMIT=$time_lim"
                cmd="$cmd bash submit_asso_real.sh"
            else
                cmd="NSEEDS=$NSEEDS TIMEOUT_OVERRIDE=$budget REGIME=raw"
                cmd="$cmd RUN_NAME=$run_name TIME_LIMIT=$time_lim"
                cmd="$cmd bash submit_asso_table2.sh"
            fi
            ;;
    esac

    echo "=== [$run_name] time=$time_lim ==="
    echo "$cmd"
    if [ $DRY_RUN -eq 0 ]; then
        eval "$cmd" || echo "[warn] $run_name a echoue"
    fi
    echo ""
}

for method in $METHODS; do
    [ -n "$METHOD_FILTER" ] && [ "$method" != "$METHOD_FILTER" ] && continue
    for rank in $RANKS; do
        [ -n "$RANK_FILTER" ] && [ "$rank" != "$RANK_FILTER" ] && continue
        for budget in $BUDGETS; do
            [ -n "$BUDGET_FILTER" ] && [ "$budget" != "$BUDGET_FILTER" ] && continue
            launch "$method" "$rank" "$budget"
        done
    done
done

echo ""
echo ">>> Fait. Suivre : squeue -u \$USER"
echo ">>> Repertoires : cluster/results_<method>_<rank>_<budget>s/"
