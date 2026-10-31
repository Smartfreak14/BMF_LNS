#!/bin/bash
# Etat des campagnes de launch_rerun.sh : par campagne, taches terminees / attendues
# (ligne "Finished" dans logs/<campagne>_<job>_<tache>.out) et jobs encore dans la file.
# Une campagne sans aucun log est "non soumise", ou "soumise" (submitted_jobs.log) avec son etat sacct.
# Quand plus rien n'est en file et que toutes les taches sont terminees, genere le tableau
# de resultats table_results.tex (make_frombelow_table.py).
#
# Usage (depuis cluster/) :  bash status_rerun.sh          (ONLY="..." pour restreindre)
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

CAMPAIGNS="boolmaxls_raw_k:datasets_real.csv boolmaxls_raw_khalf:datasets_real_khalf.csv
boolmaxls_simpluni:datasets.csv boolmaxls_simplexi:datasets_simplexi.csv
ablation_raw_k:datasets_real.csv ablation_raw_khalf:datasets_real_khalf.csv
nuwls_greedy_raw_k:datasets_real.csv nuwls_greedy_raw_khalf:datasets_real_khalf.csv
boolmaxls_simpluni_oldk:datasets_simpluni_oldk.csv boolmaxls_simplexi_oldk:datasets_simplexi_oldk.csv
nuwls_greedy_raw_k_large:datasets_real_nuwls_large.csv nuwls_greedy_raw_khalf_large:datasets_real_khalf_nuwls_large.csv"
for a in grecond greess mdlgrecond; do
    CAMPAIGNS+=" ${a}_raw_k:datasets_real.csv ${a}_raw_khalf:datasets_real_khalf.csv"
    CAMPAIGNS+=" ${a}_simpluni:datasets.csv ${a}_simplexi:datasets_simplexi.csv ${a}_synthetic:datasets_synthetic.csv"
done
# Configuration finale : <campagne>:<liste>:<taches par matrice> (une par graine).
S="${NSEEDS_FINAL:-3}"
for c in boolmaxls_raw_k_final:datasets_real.csv boolmaxls_raw_khalf_final:datasets_real_khalf.csv \
         nuwls_greedy_raw_k_final:datasets_real.csv nuwls_greedy_raw_khalf_final:datasets_real_khalf.csv \
         boolmaxls_synthetic:datasets_synthetic.csv nuwls_greedy_synthetic:datasets_synthetic.csv \
         boolmaxls_greess_raw_k:datasets_real.csv boolmaxls_greess_raw_khalf:datasets_real_khalf.csv \
         boolmaxls_nodf_raw_k:datasets_real.csv boolmaxls_nodf_raw_khalf:datasets_real_khalf.csv \
         boolmaxls_raw_k_final_redo:datasets_real_redo.csv boolmaxls_raw_khalf_final_redo:datasets_real_khalf_redo.csv; do
    CAMPAIGNS+=" $c:$S"
done

QUEUED=0
if command -v squeue > /dev/null; then
    QUEUED=$(squeue -u "$USER" -h -t PENDING,RUNNING,CONFIGURING,COMPLETING | wc -l)
    echo "Jobs dans la file : $QUEUED  ($(squeue -u "$USER" -h -t RUNNING | wc -l) en cours)"
fi

ALL_DONE=1
printf "%-30s %9s  %s\n" campagne terminees etat
for c in $CAMPAIGNS; do
    run="${c%%:*}"; list="${c#*:}"; per=1
    [[ "$list" == *:* ]] && { per="${list#*:}"; list="${list%%:*}"; }
    if [ -n "${ONLY:-}" ] && [[ " $ONLY " != *" $run "* ]]; then continue; fi
    n=$(( $(awk 'NR==1 {next} /^#/ {next} NF {x++} END {print x+0}' "$list") * per ))
    # <campagne>_<job>_<tache>.out : le chiffre apres le nom exclut <campagne>_final, _large, etc.
    done=$(grep -l "^Finished" logs/"${run}"_[0-9]*_[0-9]*.out 2> /dev/null | sed -E 's/.*_([0-9]+)\.out$/\1/' | sort -u | wc -l)
    if ! ls logs/"${run}"_[0-9]*_[0-9]*.out > /dev/null 2>&1; then
        # Aucune tache demarree : jamais soumise, ou soumise et encore en file ou rejetee.
        id=$(awk -v r="$run" '$1 == r {id = $2} END {print id}' submitted_jobs.log 2> /dev/null)
        if [ -z "$id" ]; then
            state="non soumise"
        else
            st=$(sacct -j "$id" -X -n -P -o State 2> /dev/null | awk '{print $1}' | sort -u | paste -sd, -)
            state="soumise (job $id), aucune tache demarree : ${st:-etat inconnu de sacct}"
        fi
        ALL_DONE=0
    elif [ "$done" -ge "$n" ]; then
        state="finie"
    elif [ "$QUEUED" -gt 0 ]; then
        state="en cours"; ALL_DONE=0
    else
        missing=$(comm -23 <(seq 1 "$n" | sort) <(grep -l "^Finished" logs/"${run}"_[0-9]*_[0-9]*.out 2> /dev/null \
                  | sed -E 's/.*_([0-9]+)\.out$/\1/' | sort -u) | sort -n | tr '\n' ' ')
        state="INCOMPLETE : taches $missing(voir logs/${run}_*.err, relancer avec ONLY=$run)"; ALL_DONE=0
    fi
    printf "%-30s %4d/%-4d  %s\n" "$run" "$done" "$n" "$state"
done

if [ "$ALL_DONE" = 1 ] && [ "$QUEUED" = 0 ]; then
    echo ""
    echo ">>> Toutes les campagnes sont finies."
    python3 make_frombelow_table.py -o table_results.tex && echo ">>> Tableau : $HERE/table_results.tex"
else
    echo ""
    echo ">>> Pas encore fini. Relancez : bash status_rerun.sh"
fi
