#!/bin/bash
# Relance complete BoolMaxLS (budget uniforme 30 min, NSEEDS graines) :
#   boolmaxls_raw_k       : 31 matrices brutes, k        (datasets_real.csv)
#   boolmaxls_raw_khalf   : 31 matrices brutes, floor(k/2) (datasets_real_khalf.csv)
#   boolmaxls_simpluni    : Simpli∀  (Real/SimplUni, datasets.csv)
#   boolmaxls_simplexi    : Simpli∃₀ (Real/SimplExi, datasets_simplexi.csv)
#   ablation_raw_k        : 5 variantes sur les brutes a k
#   ablation_raw_khalf    : 5 variantes sur les brutes a floor(k/2)
#   nuwls_greedy_raw_k     : NuWLS-c demarre du glouton, brutes a k        (graine 1)
#   nuwls_greedy_raw_khalf : NuWLS-c demarre du glouton, brutes a floor(k/2) (graine 1)
#   <algo>_raw_k / <algo>_raw_khalf   : baselines from-below deterministes (algo = grecond,
#       greess, mdlgrecond) sur les brutes, k ou floor(k/2) facteurs au plus
#   <algo>_simpluni / <algo>_simplexi : memes baselines jusqu'a la factorisation exacte
#   (cellules inconnues : free, puis zero et one si la matrice en a)
# Complements, lances seulement par ONLY :
#   boolmaxls_simpluni_oldk / boolmaxls_simplexi_oldk : Simpli aux rangs de l'ancien tableau
#       absents de datasets.csv / datasets_simplexi.csv (datasets_simpl*_oldk.csv)
#   nuwls_greedy_raw_k_large / nuwls_greedy_raw_khalf_large : NuWLS-c sur les encodages au-dessus
#       de MAX_CLAUSES, jusqu'a MAX_CLAUSES_LARGE (datasets_real*_nuwls_large.csv)
# Configuration finale, NSEEDS_FINAL graines a raison d'une tache SLURM par (matrice, graine),
# lancee seulement par ONLY :
#   boolmaxls_raw_k_final / boolmaxls_raw_khalf_final : BoolMaxLS sur les brutes, k et floor(k/2)
#   nuwls_greedy_raw_k_final / nuwls_greedy_raw_khalf_final : NuWLS-c demarre du glouton, limite
#       MAX_CLAUSES_LARGE (remplace les campagnes _large)
#   boolmaxls_synthetic / nuwls_greedy_synthetic : les 12 instances synthetiques (datasets_synthetic.csv)
#   boolmaxls_greess_raw_k / boolmaxls_greess_raw_khalf : BoolMaxLS demarre des facteurs de GreEss
#   <algo>_synthetic : baselines from-below sur les 12 instances synthetiques, au rang k plante
#   boolmaxls_nodf_raw_k / boolmaxls_nodf_raw_khalf : BoolMaxLS sans double flip (--ablation
#       no-double-flip), memes reglages que boolmaxls_raw_k_final / _khalf_final
#   boolmaxls_raw_k_final_redo / boolmaxls_raw_khalf_final_redo : relance des matrices de
#       datasets_real*_redo.csv, dont les executions _final ont ete tuees ; leurs resultats
#       remplacent ceux de _final dans les tableaux
#   Ces quatre campagnes coupent l'appel LNS en cours a la fin du budget (--lns-time-limit inf) :
#   meme erreur a 30 min que sans limite, mais chaque execution finit et est verifiee.
#
# TIME_LIMIT = NSEEDS x (runs chronometres x BUDGET + runs x VERIF_MARGIN) + START_MARGIN.
# Configuration finale : TL_FINAL = 6 x BUDGET par graine pour BoolMaxLS, dont le sous-probleme exact
# du LNS n'a pas de limite de temps (5 taches coupees a 50 min dans la relance) ; NuWLS-c, arrete
# au budget par bench_maxsat.sh, et BoolMaxLS avec --lns-time-limit gardent
# BUDGET + VERIF_MARGIN + START_MARGIN.
# Les autres variables lues par les jobs (racine, donnees, binaire NuWLS, traitements des
# inconnues) sont fixees ici : une variable exportee dans la session ne les change pas.
# Une campagne dont le job est encore en file (submitted_jobs.log, squeue) n'est pas resoumise.
#
# Env : NSEEDS (1), NSEEDS_FINAL (3), BUDGET (1800), VERIF_MARGIN (300), START_MARGIN (900),
#       MAX_CLAUSES (20000000, encodages plus gros non lances), MAX_CLAUSES_LARGE (60000000),
#       ONLY (liste de campagnes, defaut toutes sauf les complements), DRY_RUN=1 (affiche sans soumettre).
#
# Usage (noeud de login, depuis cluster/) :  bash launch_rerun.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
CPP_ROOT="$(cd "$HERE/.." && pwd)"

NSEEDS="${NSEEDS:-1}"
NSEEDS_FINAL="${NSEEDS_FINAL:-3}"
BUDGET="${BUDGET:-1800}"
VERIF_MARGIN="${VERIF_MARGIN:-300}"
START_MARGIN="${START_MARGIN:-900}"
MAX_CLAUSES="${MAX_CLAUSES:-20000000}"
MAX_CLAUSES_LARGE="${MAX_CLAUSES_LARGE:-60000000}"
FROMBELOW_ALL=""
for a in grecond greess mdlgrecond; do
    FROMBELOW_ALL+=" ${a}_raw_k ${a}_raw_khalf ${a}_simpluni ${a}_simplexi"
done
ONLY="${ONLY:-boolmaxls_raw_k boolmaxls_raw_khalf boolmaxls_simpluni boolmaxls_simplexi ablation_raw_k ablation_raw_khalf nuwls_greedy_raw_k nuwls_greedy_raw_khalf $FROMBELOW_ALL}"
DRY_RUN="${DRY_RUN:-0}"

hms() { local m=$(( ($1 + 59) / 60 )); printf "%02d:%02d:00" $(( m / 60 )) $(( m % 60 )); }

# 1 run chronometre par graine ; ablation : 4 chronometres (greedy instantane) + 5 verifs.
TL_BOOLMAXLS=$(hms $(( NSEEDS * (BUDGET + VERIF_MARGIN) + START_MARGIN )))
TL_ABLATION=$(hms $(( NSEEDS * (4 * BUDGET + 5 * VERIF_MARGIN) + START_MARGIN )))
# Baselines from-below : 3 traitements des inconnues x (1 calcul + 1 verif).
TL_FROMBELOW=$(hms $(( 3 * (BUDGET + VERIF_MARGIN) + START_MARGIN )))
TL_FINAL=$(hms $(( 6 * BUDGET )))
TL_NUWLS_SEED=$(hms $(( BUDGET + VERIF_MARGIN + START_MARGIN )))
TL_BOUNDED="$TL_NUWLS_SEED"

[ -x "$CPP_ROOT/build/BoolMaxLS" ] || { echo "ERREUR : build/BoolMaxLS absent (bash setup_env.sh)"; exit 1; }
[ -f "$CPP_ROOT/src/verif.py" ]    || { echo "ERREUR : src/verif.py absent"; exit 1; }

# Une ligne par campagne soumise (campagne, job, date), ecrite par les submit_*.sh et lue
# par status_rerun.sh.
SUBMITTED="$HERE/submitted_jobs.log"
touch "$SUBMITTED"
NSUB0=$(wc -l < "$SUBMITTED")

boolmaxls() {
    local run="$1" list="$2" base="$3" regime="$4"
    echo ">>> $run : $list  data/$base  regime=$regime  budget=${BUDGET}s  seeds=$NSEEDS  time=$TL_BOOLMAXLS"
    [ "$DRY_RUN" = 1 ] && return
    RUN_NAME="$run" DATASETS_FILE="$HERE/$list" DATA_BASE="$base" REGIME="$regime" \
    NSEEDS="$NSEEDS" SEED_PER_TASK=0 INIT_ALGO="" BOOLMAXLS_ARGS="" TIMEOUT_OVERRIDE="$BUDGET" TIME_LIMIT="$TL_BOOLMAXLS" \
        bash "$HERE/submit_boolmaxls_real.sh"
}

# Configuration finale : une tache par (matrice, graine) ; init = algorithme de depart, sinon le glouton ;
# args = options de BoolMaxLS. Avec --lns-time-limit, chaque execution s'arrete au budget : TL_BOUNDED.
boolmaxls_final() {
    local run="$1" list="$2" base="$3" regime="$4" init="${5:-}" args="${6:-}" tl="$TL_FINAL"
    if [[ " $args " == *" --lns-time-limit "* ]]; then
        tl="$TL_BOUNDED"
        [[ "$("$CPP_ROOT/build/BoolMaxLS" --help)" == *--lns-time-limit* ]] \
            || { echo "ERREUR : build/BoolMaxLS ne connait pas --lns-time-limit (bash setup_env.sh)"; exit 1; }
    fi
    echo ">>> $run : $list  data/$base  regime=$regime  depart=${init:-glouton}${args:+  options=$args}  budget=${BUDGET}s  graines=$NSEEDS_FINAL (une tache chacune)  time=$tl"
    [ "$DRY_RUN" = 1 ] && return
    RUN_NAME="$run" DATASETS_FILE="$HERE/$list" DATA_BASE="$base" REGIME="$regime" \
    NSEEDS="$NSEEDS_FINAL" SEED_PER_TASK=1 INIT_ALGO="$init" BOOLMAXLS_ARGS="$args" TIMEOUT_OVERRIDE="$BUDGET" TIME_LIMIT="$tl" \
        bash "$HERE/submit_boolmaxls_real.sh"
}

ablation() {
    local run="$1" list="$2"
    echo ">>> $run : $list  data/Real  5 variantes  budget=${BUDGET}s  seeds=$NSEEDS  time=$TL_ABLATION"
    [ "$DRY_RUN" = 1 ] && return
    (
        export PROJECT_ROOT="$CPP_ROOT" JOB_SCRIPT="job_ablation_array.sh" NSEEDS SEED_PER_TASK=0 TIMEOUT_OVERRIDE="$BUDGET" \
               REGIME="raw" TIME_LIMIT="$TL_ABLATION" MEM_PER_TASK="${MEM_ABLATION:-16G}"
        source "$HERE/_submit_common.sh"
        submit_array "$HERE/$list" "Real" "$run"
    )
}

# nseeds vide : graine 1 seulement ; sinon une tache par (matrice, graine).
nuwls() {
    local run="$1" list="$2" max="$3" nseeds="${4:-}" base="${5:-Real}"
    local spt=0 tl="$TL_BOOLMAXLS" seeds="seed=1"
    [ -z "$nseeds" ] || { spt=1; tl="$TL_NUWLS_SEED"; seeds="graines=$nseeds (une tache chacune)"; }
    echo ">>> $run : $list  data/$base  NuWLS-c depart glouton  budget=${BUDGET}s  $seeds  max=${max} clauses  time=$tl"
    [ "$DRY_RUN" = 1 ] && return
    [ -x "$CPP_ROOT/third_party/NuWLS/nuwls_c" ] || { echo "ERREUR : third_party/NuWLS/nuwls_c absent"; exit 1; }
    (
        export PROJECT_ROOT="$CPP_ROOT" JOB_SCRIPT="job_maxsat_array.sh" MAXSAT_TIMEOUT="$BUDGET" \
               NUWLS_INIT="greedy" NUWLS_BIN="$CPP_ROOT/third_party/NuWLS/nuwls_c" SEED=1 MAX_CLAUSES="$max" \
               SEED_PER_TASK="$spt" NSEEDS="${nseeds:-1}" TIME_LIMIT="$tl" MEM_PER_TASK="${MEM_NUWLS:-32G}"
        source "$HERE/_submit_common.sh"
        submit_array "$HERE/$list" "$base" "$run"
    )
}

frombelow() {
    local run="$1" list="$2" base="$3" exact="$4"
    local algo="${run%%_*}"
    case "$algo" in grecond|greess|mdlgrecond) ;; *) echo "ERREUR : campagne inconnue $run"; exit 1 ;; esac
    echo ">>> $run : $list  data/$base  $algo  $([ "$exact" = 1 ] && echo "jusqu'a l'exact" || echo "k de la liste")  budget=${BUDGET}s  inconnues=free,zero,one  time=$TL_FROMBELOW"
    [ "$DRY_RUN" = 1 ] && return
    [ -x "$CPP_ROOT/build/FromBelowBaselines" ] || { echo "ERREUR : build/FromBelowBaselines absent (bash setup_env.sh)"; exit 1; }
    (
        export PROJECT_ROOT="$CPP_ROOT" JOB_SCRIPT="job_frombelow_dataset.sh" ALGO="$algo" EXACT="$exact" \
               MISSING_MODES="free zero one" BUDGET SEED_PER_TASK=0 TIME_LIMIT="$TL_FROMBELOW" \
               MEM_PER_TASK="${MEM_FROMBELOW:-4G}"
        source "$HERE/_submit_common.sh"
        submit_array "$HERE/$list" "$base" "$run"
    )
}

for c in $ONLY; do
    # Deux jobs de la meme campagne ecriraient dans le meme results_<campagne>/.
    id=$(awk -v r="$c" '$1 == r {id = $2} END {print id}' "$SUBMITTED")
    if [ -n "$id" ] && [ -n "$(squeue -h -j "$id" 2> /dev/null)" ]; then
        echo ">>> $c : deja en file (job $id), non resoumise"
        continue
    fi
    case "$c" in
        boolmaxls_raw_k)     boolmaxls "$c" datasets_real.csv       Real          raw ;;
        boolmaxls_raw_khalf) boolmaxls "$c" datasets_real_khalf.csv Real          raw ;;
        boolmaxls_simpluni)  boolmaxls "$c" datasets.csv            Real/SimplUni simpluni ;;
        boolmaxls_simplexi)  boolmaxls "$c" datasets_simplexi.csv   Real/SimplExi simplexi ;;
        boolmaxls_simpluni_oldk) boolmaxls "$c" datasets_simpluni_oldk.csv Real/SimplUni simpluni ;;
        boolmaxls_simplexi_oldk) boolmaxls "$c" datasets_simplexi_oldk.csv Real/SimplExi simplexi ;;
        ablation_raw_k)      ablation "$c" datasets_real.csv ;;
        ablation_raw_khalf)  ablation "$c" datasets_real_khalf.csv ;;
        nuwls_greedy_raw_k)     nuwls "$c" datasets_real.csv       "$MAX_CLAUSES" ;;
        nuwls_greedy_raw_khalf) nuwls "$c" datasets_real_khalf.csv "$MAX_CLAUSES" ;;
        nuwls_greedy_raw_k_large)     nuwls "$c" datasets_real_nuwls_large.csv       "$MAX_CLAUSES_LARGE" ;;
        nuwls_greedy_raw_khalf_large) nuwls "$c" datasets_real_khalf_nuwls_large.csv "$MAX_CLAUSES_LARGE" ;;
        boolmaxls_raw_k_final)        boolmaxls_final "$c" datasets_real.csv       Real      raw ;;
        boolmaxls_raw_khalf_final)    boolmaxls_final "$c" datasets_real_khalf.csv Real      raw ;;
        boolmaxls_synthetic)          boolmaxls_final "$c" datasets_synthetic.csv  Synthetic synth ;;
        boolmaxls_greess_raw_k)       boolmaxls_final "$c" datasets_real.csv       Real      raw greess ;;
        boolmaxls_greess_raw_khalf)   boolmaxls_final "$c" datasets_real_khalf.csv Real      raw greess ;;
        boolmaxls_nodf_raw_k)         boolmaxls_final "$c" datasets_real.csv       Real      raw "" "--ablation no-double-flip --lns-time-limit inf" ;;
        boolmaxls_nodf_raw_khalf)     boolmaxls_final "$c" datasets_real_khalf.csv Real      raw "" "--ablation no-double-flip --lns-time-limit inf" ;;
        boolmaxls_raw_k_final_redo)     boolmaxls_final "$c" datasets_real_redo.csv       Real raw "" "--lns-time-limit inf" ;;
        boolmaxls_raw_khalf_final_redo) boolmaxls_final "$c" datasets_real_khalf_redo.csv Real raw "" "--lns-time-limit inf" ;;
        nuwls_greedy_raw_k_final)     nuwls "$c" datasets_real.csv       "$MAX_CLAUSES_LARGE" "$NSEEDS_FINAL" ;;
        nuwls_greedy_raw_khalf_final) nuwls "$c" datasets_real_khalf.csv "$MAX_CLAUSES_LARGE" "$NSEEDS_FINAL" ;;
        nuwls_greedy_synthetic)       nuwls "$c" datasets_synthetic.csv  "$MAX_CLAUSES_LARGE" "$NSEEDS_FINAL" Synthetic ;;
        *_raw_k)        frombelow "$c" datasets_real.csv       Real          0 ;;
        *_raw_khalf)    frombelow "$c" datasets_real_khalf.csv Real          0 ;;
        *_simpluni)     frombelow "$c" datasets.csv            Real/SimplUni 1 ;;
        *_simplexi)     frombelow "$c" datasets_simplexi.csv   Real/SimplExi 1 ;;
        *_synthetic)    frombelow "$c" datasets_synthetic.csv  Synthetic     0 ;;
        *) echo "ERREUR : campagne inconnue $c"; exit 1 ;;
    esac
done

echo ""
if [ "$DRY_RUN" = 1 ]; then
    echo ">>> DRY_RUN=1 : simulation, aucune campagne soumise. Pour lancer : bash launch_rerun.sh"
else
    echo ">>> $(( $(wc -l < "$SUBMITTED") - NSUB0 )) campagne(s) soumise(s) (voir submitted_jobs.log)."
    echo ">>> Suivi : sq ; bash status_rerun.sh"
fi
