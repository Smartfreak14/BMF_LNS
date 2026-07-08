#!/bin/bash
# Compile les solveurs baseline MaxSAT anytime dans third_party/.
# EvalMaxSAT2022 est deja compile via setup_env.sh.
#
# Cible principale : cluster Narval (Linux + gcc + libs statiques).
# Sur Mac, `-static` echoue souvent : on l'enleve, mais des solveurs peuvent
# quand meme rater. Pour un vrai benchmark, faire tourner ce script sur le
# cluster, pas sur Mac.

set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
CPP_ROOT="$(cd "$HERE/.." && pwd)"
THIRD_PARTY="$CPP_ROOT/third_party"
mkdir -p "$THIRD_PARTY"

IS_MAC=0
[ "$(uname)" = "Darwin" ] && IS_MAC=1

try_clone() {
    local name="$1"; shift
    local dir="$THIRD_PARTY/$name"
    [ -d "$dir" ] && { echo "  [$name] deja present, skip clone"; return 0; }
    for url in "$@"; do
        echo "  [$name] git clone $url ..."
        if git clone --depth 1 "$url" "$dir" 2>/dev/null; then
            echo "  [$name] clone OK"
            return 0
        fi
    done
    echo "  [$name] echec du clone (URLs testees : $*)"
    echo "  [$name] telecharger manuellement, poser dans third_party/$name/, relancer"
    return 1
}

# =====================================================================
echo ">>> EvalMaxSAT (Avellaneda, MaxSAT Eval champion complet)"
if [ ! -x "$THIRD_PARTY/EvalMaxSAT/build/main/EvalMaxSAT_bin" ]; then
    if [ ! -d "$THIRD_PARTY/EvalMaxSAT" ]; then
        (cd "$THIRD_PARTY" && git clone --depth 1 \
            https://github.com/FlorentAvellaneda/EvalMaxSAT.git)
    fi
    (cd "$THIRD_PARTY/EvalMaxSAT" && mkdir -p build && cd build && cmake .. && make -j4)
    if [ -x "$THIRD_PARTY/EvalMaxSAT/build/main/EvalMaxSAT_bin" ]; then
        echo "  [EvalMaxSAT] OK"
    else
        echo "  [EvalMaxSAT] echec du build"
    fi
else
    echo "  [EvalMaxSAT] deja compile, skip"
fi
echo ""

# =====================================================================
echo ">>> NuWLS (Chu, Cai, Luo, AAAI 2023)"
if try_clone NuWLS "https://github.com/filyouzicha/NuWLS.git"; then
    # NuWLS a une sous-arbo NuWLS-source-code/ avec un Makefile
    NUWLS_SRC="$THIRD_PARTY/NuWLS/NuWLS-source-code"
    if [ -d "$NUWLS_SRC" ] && [ -f "$NUWLS_SRC/Makefile" ]; then
        # Sur Mac : virer -static (crt0.o absent) ; sur Linux : garder.
        if [ $IS_MAC -eq 1 ]; then
            sed -i.bak 's/-static//' "$NUWLS_SRC/Makefile"
        fi
        (cd "$NUWLS_SRC" && make -j4) && \
        # Le binaire s'appelle 'nuwls', on symlink en 'nuwls_c' pour bench_maxsat.sh
        [ -x "$NUWLS_SRC/nuwls" ] && \
            ln -sf "$NUWLS_SRC/nuwls" "$THIRD_PARTY/NuWLS/nuwls_c" && \
            echo "  [NuWLS] OK -> third_party/NuWLS/nuwls_c"
    else
        echo "  [NuWLS] structure inattendue, Makefile introuvable"
    fi
fi
echo ""

# =====================================================================
echo ">>> SATLike-c (Lei, Cai, MaxSAT Eval 2021)"
# Pas d'URL canonique connue ; le code est distribue sur la page MaxSAT
# Evaluation 2021. Pour l'installer manuellement :
if [ ! -x "$THIRD_PARTY/SATLike-c/satlike_c" ]; then
    cat <<'EOF'
  [SATLike-c] installation manuelle requise :
    1. Telecharger la soumission SATLike-c depuis MaxSAT Evaluation 2021
       https://maxsat-evaluations.github.io/2021/mse21-solver-src/incomplete/SATLike-c.zip
       (URL a verifier, peut avoir change)
    2. Extraire dans third_party/SATLike-c/
    3. cd third_party/SATLike-c && make
    4. Le binaire attendu est satlike_c
EOF
fi
echo ""

# =====================================================================
echo ">>> loandra (source deja presente)"
if [ -d "$THIRD_PARTY/loandra" ] && [ -f "$THIRD_PARTY/loandra/Makefile" ]; then
    if [ $IS_MAC -eq 1 ]; then
        echo "  [loandra] SKIP sur Mac (architecture mismatch gmp x86_64/arm64)"
        echo "           Fonctionnera sur Linux/cluster."
    else
        (cd "$THIRD_PARTY/loandra" && make -j4 r)
    fi
else
    echo "  [loandra] source absente"
fi
echo ""

# =====================================================================
echo ""
echo ">>> Bilan des binaires disponibles :"
for entry in \
    "EvalMaxSAT/build/main/EvalMaxSAT_bin" \
    "NuWLS/nuwls_c" \
    "SATLike-c/satlike_c" \
    "loandra/loandra"; do
    full="$THIRD_PARTY/$entry"
    if [ -x "$full" ]; then
        echo "   OK  $entry"
    else
        echo "   NON $entry"
    fi
done

if [ $IS_MAC -eq 1 ]; then
    echo ""
    echo ">>> Note MacOS : les compilations foireuses sont attendues (-static, gmp)."
    echo ">>> Pour le benchmark AAAI, executer ce script sur Narval (Linux)."
fi
