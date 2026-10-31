#!/bin/bash
# ============================================================
#  setup_env.sh : compilation complete du projet BoolMaxLS
#  sur Calcul Québec / Alliance Canada (cluster Narval).
#
#  - Refuse de compiler si les sources ne sont pas celles copiees par
#    deploy_narval.sh (SOURCES.sha256), ou si des campagnes de ce dossier
#    sont encore en file (submitted_jobs.log, squeue ; FORCE=1 pour passer outre).
#  - Charge les modules StdEnv / gcc / cmake.
#  - Compile EvalMaxSAT2022 (third_party/, makefile interne) si la
#    bibliotheque manque ou n'est pas un objet Linux (ELF).
#  - Recompile de zero BoolMaxLS et FromBelowBaselines (CMake, Release)
#    et NuWLS-c (third_party/NuWLS/nuwls_c).
#
#  Tout le code et les binaires restent dans /project/def-florenta
#  (quota partagé). Pas de venv Python ici : aggregate.py n'utilise
#  que la stdlib donc le 'python' fourni par le module suffit.
#
#  À lancer DEPUIS le nœud de login :
#      cd <projet>/cluster
#      bash setup_env.sh
#
#  Compile toujours le projet qui contient ce script (dossier parent de
#  cluster/) : un PROJECT_ROOT exporte dans la session est ignore.
# ============================================================
set -euo pipefail

CPP_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if [ -n "${PROJECT_ROOT:-}" ] && [ "$(cd "$PROJECT_ROOT" 2>/dev/null && pwd)" != "$CPP_ROOT" ]; then
    echo ">>> ATTENTION : PROJECT_ROOT=$PROJECT_ROOT (session) ignore."
fi
[ -f "$CPP_ROOT/CMakeLists.txt" ] || { echo "ERREUR : pas de CMakeLists.txt sous $CPP_ROOT"; exit 1; }
EVALMAXSAT_DIR="$CPP_ROOT/third_party/EvalMaxSAT2022"
BUILD_DIR="$CPP_ROOT/build"

echo ">>> Projet : $CPP_ROOT"
echo ""

# Sources identiques a celles du Mac (SOURCES.sha256, ecrit par deploy_narval.sh),
# verifiees avant toute compilation et avant d'effacer build/.
if [ ! -f "$CPP_ROOT/SOURCES.sha256" ]; then
    echo "ERREUR : pas de SOURCES.sha256 dans $CPP_ROOT : ce dossier n'a pas ete copie"
    echo "         par deploy_narval.sh (sur le Mac : bash cluster/deploy_narval.sh)."
    exit 1
fi
if command -v sha256sum >/dev/null; then SHA=(sha256sum); else SHA=(shasum -a 256); fi
( cd "$CPP_ROOT" && "${SHA[@]}" --quiet -c SOURCES.sha256 ) \
    || { echo "ERREUR : fichiers differents de la copie du Mac (ci-dessus) : relancer deploy_narval.sh"; exit 1; }
echo ">>> Sources : $(wc -l < "$CPP_ROOT/SOURCES.sha256") fichiers identiques a la copie du Mac"
echo ""

# Les taches qui demarreraient pendant la recompilation ne trouveraient pas les binaires.
SUBMITTED="$CPP_ROOT/cluster/submitted_jobs.log"
if [ "${FORCE:-0}" != 1 ] && [ -f "$SUBMITTED" ] && command -v squeue > /dev/null; then
    busy=""
    for id in $(awk '{print $2}' "$SUBMITTED" | sort -u); do
        if [ -n "$(squeue -h -j "$id" 2> /dev/null)" ]; then busy+=" $id"; fi
    done
    if [ -n "$busy" ]; then
        echo "ERREUR : jobs de ce dossier encore en file :$busy"
        echo "         Les binaires sont deja compiles ; recompiler maintenant ferait echouer les"
        echo "         taches qui demarrent. Attendre leur fin, ou FORCE=1 bash setup_env.sh."
        exit 1
    fi
fi

echo ">>> Vérification quotas (ignorer si la commande n'existe pas)..."
diskusage_report 2>/dev/null || true
echo ""

echo ">>> Chargement des modules Narval..."
module --force purge
module load StdEnv/2023
module load gcc/12.3
module load cmake/3.27.7
module load python/3.11   # juste pour aggregate.py + verif.py
echo ""

echo ">>> Modules charges :"
module list 2>&1 | sed 's/^/    /'
echo ""

# --------------------------------------------------------------
# 1) Compilation d'EvalMaxSAT2022 (libipamirEvalMaxSAT2022.a)
# --------------------------------------------------------------
if [ ! -d "$EVALMAXSAT_DIR" ]; then
    echo "ERREUR : EvalMaxSAT2022 manquant dans $EVALMAXSAT_DIR"
    echo "         Il doit etre commit avec le projet (third_party/)."
    exit 1
fi

LIB="$EVALMAXSAT_DIR/libipamirEvalMaxSAT2022.a"

# La lib peut etre presente mais venir d'un upload depuis un Mac (Mach-O) :
# ld Linux ne resoudrait pas les symboles. On verifie qu'elle contient bien
# un objet ELF avant de l'accepter.
need_rebuild=1
if [ -f "$LIB" ]; then
    if ar t "$LIB" >/dev/null 2>&1 && \
       ar x --output=/tmp "$LIB" "$(ar t "$LIB" | head -1)" 2>/dev/null && \
       file -b "/tmp/$(ar t "$LIB" | head -1)" 2>/dev/null | grep -q '^ELF'; then
        need_rebuild=0
        echo ">>> EvalMaxSAT2022 deja compile (ELF) : $LIB"
    else
        echo ">>> $LIB present mais pas ELF (build d'une autre plateforme) -> rebuild."
        rm -f "$LIB"
    fi
    rm -f "/tmp/$(ar t "$LIB" 2>/dev/null | head -1)" 2>/dev/null || true
fi

if [ "$need_rebuild" -eq 1 ]; then
    echo ">>> Compilation EvalMaxSAT2022 (peut prendre 5-15 min)..."
    (
        cd "$EVALMAXSAT_DIR"
        make clean 2>/dev/null || true
        find . -name '*.o' -delete
        make -j"${SLURM_CPUS_ON_NODE:-4}"
    )
    if [ ! -f "$LIB" ]; then
        echo "ERREUR FATALE : libipamirEvalMaxSAT2022.a non produit."
        exit 1
    fi
    echo "    EvalMaxSAT2022 compile : $LIB"
fi
echo ""

# --------------------------------------------------------------
# 2) Configuration + compilation BoolMaxLS
# --------------------------------------------------------------
# Build de zero : aucun objet d'une compilation precedente n'est reutilise.
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"
(
    cd "$BUILD_DIR"
    echo ">>> CMake configure..."
    cmake .. -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_C_COMPILER=gcc \
             -DCMAKE_CXX_COMPILER=g++
    echo ""
    echo ">>> CMake build..."
    make -j"${SLURM_CPUS_ON_NODE:-4}"
)

# --------------------------------------------------------------
# 3) Compilation NuWLS-c (nuwls_c, utilise par scripts/bench_maxsat.sh)
# --------------------------------------------------------------
NUWLS_DIR="$CPP_ROOT/third_party/NuWLS/NuWLS-source-code"
echo ""
echo ">>> Compilation NuWLS-c..."
( cd "$NUWLS_DIR" && rm -f nuwls && make )
ln -sfn "$NUWLS_DIR/nuwls" "$CPP_ROOT/third_party/NuWLS/nuwls_c"

BIN="$BUILD_DIR/BoolMaxLS"
FB_BIN="$BUILD_DIR/FromBelowBaselines"
NUWLS_BIN="$CPP_ROOT/third_party/NuWLS/nuwls_c"
echo ""
echo ">>> Binaires :"
for b in "$BIN" "$FB_BIN" "$NUWLS_BIN"; do
    [ -x "$b" ] || { echo "ERREUR FATALE : $b non produit."; exit 1; }
    echo "    $(basename "$b") : $(file -b -L "$b" | cut -d, -f1-2)"
done

echo ""
echo ">>> Verification : ./BoolMaxLS --help"
"$BIN" --help | head -20 | sed 's/^/    /'
echo ""

# --------------------------------------------------------------
# 4) Mini tests : BoolMaxLS + verif.py, chargeur CSV, NuWLS -init
# --------------------------------------------------------------
ZOO="$CPP_ROOT/data/Real/SimplUni/zoo.csv"
if [ -f "$ZOO" ]; then
    echo ">>> Mini test : BoolMaxLS + verif.py sur zoo.csv (k=25, timeout 5s)..."
    ( cd "$BUILD_DIR" && ./BoolMaxLS "$ZOO" 25 5 > /tmp/setup_env_zoo.log 2>&1 ) \
        && ! grep -q -E "MISMATCH|n/a" /tmp/setup_env_zoo.log \
        && echo "    OK" \
        || { echo "ERREUR : run sanity zoo a echoue, voir /tmp/setup_env_zoo.log"; exit 1; }
else
    echo ">>> Mini test skip : $ZOO introuvable."
fi

T="$(mktemp -d)"
echo ">>> Mini test : FromBelowBaselines (grecond, votes.csv, k=17)..."
"$FB_BIN" grecond "$CPP_ROOT/data/Real/votes.csv" 17 --out-dir "$T" > "$T/fb.log" 2>&1 \
    && grep -q "unknown=392" "$T/fb.log" \
    && echo "    OK (392 cellules inconnues lues)" \
    || { echo "ERREUR : FromBelowBaselines, voir $T/fb.log"; exit 1; }

echo ">>> Mini test : BoolMaxLS depart des facteurs de GreConD (votes.csv, k=17, 3s)..."
FB_ERR=$(sed -n 's/^ *k=17 factors=[0-9]* error=\([0-9]*\).*/\1/p' "$T/fb.log" | head -1)
( cd "$BUILD_DIR" && ./BoolMaxLS "$CPP_ROOT/data/Real/votes.csv" 17 3 \
    --init-factors "$T/votes_grecond_free_k17_A.csv" "$T/votes_grecond_free_k17_B.csv" > "$T/init.log" 2>&1 ) \
    && grep -q "^\[init-factors\] init_errors=$FB_ERR$" "$T/init.log" \
    && echo "    OK (depart a $FB_ERR erreurs, comme GreConD)" \
    || { echo "ERREUR : BoolMaxLS --init-factors, voir $T/init.log"; exit 1; }

# Limite 60 s : aucun appel abandonne ; limite 1 us : tous abandonnes. Dans les deux cas,
# erreur verifiee = erreur finale (colonnes 9 et 15 du CSV).
echo ">>> Mini test : BoolMaxLS --lns-time-limit (zoo.csv, k=12, 2 x 5s)..."
for LIM in 60 0.000001; do
    ( cd "$BUILD_DIR" && ./BoolMaxLS "$CPP_ROOT/data/Real/zoo.csv" 12 5 --lns-time-limit "$LIM" \
        --csv-out "$T/tl_$LIM.csv" > "$T/tl_$LIM.log" 2>&1 ) \
        && awk -F, 'NR == 2 && $9 == $15 {ok = 1} END {exit !ok}' "$T/tl_$LIM.csv" \
        || { echo "ERREUR : BoolMaxLS --lns-time-limit $LIM, voir $T/tl_$LIM.log"; exit 1; }
done
grep -q "abandonnes 0$" "$T/tl_60.log" && grep -q "abandonne apres" "$T/tl_0.000001.log" \
    && echo "    OK ($(grep '^\[LNS\]' "$T/tl_60.log"))" \
    || { echo "ERREUR : BoolMaxLS --lns-time-limit, voir $T/tl_60.log et $T/tl_0.000001.log"; exit 1; }

echo ">>> Mini test : NuWLS-c depart glouton (zoo.csv, k=12, 3s)..."
"$BIN" --dump-wcnf "$CPP_ROOT/data/Real/zoo.csv" 12 "$T/zoo.wcnf" --init-out "$T/zoo.init" > "$T/dump.log" 2>&1 \
    || { echo "ERREUR : --dump-wcnf, voir $T/dump.log"; exit 1; }
timeout -s TERM 3 "$NUWLS_BIN" "$T/zoo.wcnf" 1 -init "$T/zoo.init" > "$T/nuwls.log" 2>&1 || true
grep -q "^c init assignment" "$T/nuwls.log" && grep -q -E "^o[[:space:]]+[0-9]+" "$T/nuwls.log" \
    && echo "    OK ($(grep '^c init assignment' "$T/nuwls.log"))" \
    || { echo "ERREUR : NuWLS-c, voir $T/nuwls.log"; exit 1; }

echo ""
echo ">>> Environnement pret."
echo "    Binaires  : $BIN"
echo "                $FB_BIN"
echo "                $NUWLS_BIN"
echo "    Donnees   : $CPP_ROOT/data/Real/"
echo ""
echo "    Pour lancer la relance complete :"
echo "      cd $CPP_ROOT/cluster"
echo "      bash launch_rerun.sh"
