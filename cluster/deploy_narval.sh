#!/bin/bash
# Copie sur Narval tout ce qu'il faut pour recompiler et lancer la relance (a lancer sur le Mac).
# Aucun binaire ni objet du Mac n'est copie : setup_env.sh recompile tout sur Narval.
#   1. cree DEST sur Narval, puis y copie (une seule connexion ssh) les sources, third_party,
#      scripts, listes de jeux de donnees et donnees ;
#   2. y ajoute SOURCES.sha256 (sommes de controle locales), que setup_env.sh verifie avant de compiler.
#
# Env : REMOTE (kiko@narval.alliancecan.ca ; REMOTE= vide : copie locale, pour tester),
#       DEST (/home/kiko/projects/def-florenta/kiko/BoolMaxLS_rerun)
#
# Usage :  bash cluster/deploy_narval.sh
set -euo pipefail
CPP_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REMOTE="${REMOTE-kiko@narval.alliancecan.ca}"
DEST="${DEST:-/home/kiko/projects/def-florenta/kiko/BoolMaxLS_rerun}"

cd "$CPP_ROOT"
NUWLS=third_party/NuWLS/NuWLS-source-code
FILES=(
    CMakeLists.txt README.md src
    third_party/ipamir.h third_party/ipasir.h third_party/EvalMaxSAT2022
    $NUWLS/{pms.cpp,pms.h,basis_pms.h,build.h,deci.h,heuristic.h,Makefile}
    scripts/bench_maxsat.sh
    cluster/{_submit_common.sh,setup_env.sh,launch_rerun.sh,status_rerun.sh,deploy_narval.sh}
    cluster/{job_boolmaxls_dataset.sh,submit_boolmaxls_real.sh,job_ablation_array.sh,submit_ablation_real_khalf.sh}
    cluster/{job_maxsat_array.sh,job_frombelow_dataset.sh}
    cluster/{datasets_real.csv,datasets_real_khalf.csv,datasets.csv,datasets_simplexi.csv}
    cluster/{datasets_simpluni_oldk.csv,datasets_simplexi_oldk.csv}
    cluster/{datasets_real_nuwls_large.csv,datasets_real_khalf_nuwls_large.csv,datasets_synthetic.csv}
    cluster/{datasets_real_redo.csv,datasets_real_khalf_redo.csv}
    cluster/{analyze_boolmaxls_runs.py,analyze_ablation_khalf.py,analyze_nuwls_runs.py,analyze_frombelow_runs.py}
    cluster/{aggregate_ablation.py,aggregate_maxsat.py}
    cluster/{make_frombelow_table.py,table_main_article.tex,make_dataset_table.py,plot_traces.py}
    data/Real/*.csv data/Real/SimplUni data/Real/SimplExi data/Synthetic/*.csv
)
# Artefacts du Mac (objets, bibliotheques, en-tetes precompiles) : jamais copies.
SKIP=(.DS_Store '*.o' '*.a' '*.gch' '*.so' '*.dylib' __pycache__)

for f in "${FILES[@]}"; do
    [ -e "$f" ] || { echo "ERREUR : $CPP_ROOT/$f absent"; exit 1; }
done

EXCLUDES=()
NAMES=()
for p in "${SKIP[@]}"; do
    EXCLUDES+=(--exclude "$p")
    NAMES+=(! -name "$p")
done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
find "${FILES[@]}" -type f "${NAMES[@]}" ! -path '*/__pycache__/*' -print0 \
    | LC_ALL=C sort -z | xargs -0 shasum -a 256 > "$TMP/SOURCES.sha256"
echo ">>> $(wc -l < "$TMP/SOURCES.sha256" | tr -d ' ') fichiers -> ${REMOTE:+$REMOTE:}$DEST"

# "$TMP/./SOURCES.sha256" : avec -R, le fichier arrive a la racine de DEST.
if [ -n "$REMOTE" ]; then
    rsync -azR --rsync-path="mkdir -p '$DEST' && rsync" "${EXCLUDES[@]}" \
        "${FILES[@]}" "$TMP/./SOURCES.sha256" "$REMOTE:$DEST/"
else
    mkdir -p "$DEST"
    rsync -aR "${EXCLUDES[@]}" "${FILES[@]}" "$TMP/./SOURCES.sha256" "$DEST/"
fi

echo ">>> Copie terminee. Sur Narval :"
echo "      cd $DEST/cluster && bash setup_env.sh"
