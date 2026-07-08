#!/bin/bash
# bench_maxsat.sh : execute plusieurs solveurs MaxSAT anytime sur un WCNF.
#
# Usage :
#   bench_maxsat.sh <input.wcnf> <timeout_s> <output.csv>
#
# Chemins des solveurs (env vars, defauts pointent sur third_party/) :
#   EVALMAXSAT_BIN, NUWLS_BIN, SATLIKE_BIN, LOANDRA_BIN
# Un solveur absent est simplement skippe.

set -euo pipefail

INPUT="${1:?Usage: bench_maxsat.sh <input.wcnf> <timeout_s> <output.csv>}"
TIMEOUT="${2:?Usage: bench_maxsat.sh <input.wcnf> <timeout_s> <output.csv>}"
OUT_CSV="${3:?Usage: bench_maxsat.sh <input.wcnf> <timeout_s> <output.csv>}"

HERE="$(cd "$(dirname "$0")" && pwd)"
CPP_ROOT="$(cd "$HERE/.." && pwd)"

EVALMAXSAT_BIN="${EVALMAXSAT_BIN:-$CPP_ROOT/third_party/EvalMaxSAT/build/main/EvalMaxSAT_bin}"
NUWLS_BIN="${NUWLS_BIN:-$CPP_ROOT/third_party/NuWLS/nuwls_c}"
SATLIKE_BIN="${SATLIKE_BIN:-$CPP_ROOT/third_party/SATLike-c/satlike_c}"
LOANDRA_BIN="${LOANDRA_BIN:-$CPP_ROOT/third_party/loandra/loandra}"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

# Timeout portable (Mac/Linux) : lance en background et SIGTERM apres N secondes.
run_with_timeout() {
    local timeout="$1"; shift
    local out="$1"; shift
    "$@" > "$out" 2>&1 &
    local pid=$!
    (sleep "$timeout" && kill -TERM "$pid" 2>/dev/null && sleep 2 && kill -KILL "$pid" 2>/dev/null) &
    local killer=$!
    wait "$pid" 2>/dev/null
    local rc=$?
    kill "$killer" 2>/dev/null || true
    wait "$killer" 2>/dev/null || true
    return $rc
}

run_solver() {
    local name="$1" bin="$2"; shift 2
    # Le reste des args : gabarit d'invocation ou "@" est remplace par $INPUT.
    local args=("$@")
    if [ ! -x "$bin" ]; then
        echo "[skip] $name : $bin introuvable"
        echo "$name,N/A,0,$TIMEOUT,0,not_found" >> "$OUT_CSV"
        return
    fi

    # Substituer "@" par "$INPUT" dans les args.
    local final_args=()
    for a in "${args[@]}"; do
        if [ "$a" = "@" ]; then final_args+=("$INPUT"); else final_args+=("$a"); fi
    done

    local out="$TMPDIR/$name.out"
    local t0 t1 elapsed best solved
    t0=$(date +%s.%N)
    set +e
    run_with_timeout "$TIMEOUT" "$out" "$bin" "${final_args[@]}"
    local rc=$?
    set -e
    t1=$(date +%s.%N)
    elapsed=$(awk -v t0="$t0" -v t1="$t1" 'BEGIN{printf "%.3f", (t1-t0)*1000}')

    # Meilleur cout : derniere ligne "o <N>" (aussi "o <N> <time>" pour NuWLS).
    best=$(grep -E "^o[[:space:]]+[0-9]+" "$out" | tail -1 | awk '{print $2}' || true)
    if grep -qE "^s OPTIMUM FOUND|^c best solution found" "$out"; then solved=1
    else solved=0
    fi

    # Note explicite : crash (139=SIGSEGV, 137=SIGKILL/OOM, 143=SIGTERM/timeout).
    local note="rc=$rc"
    case "$rc" in
        139) note="crash_segv" ;;
        137) note="oom_or_killed" ;;
        143|124) note="timeout" ;;
    esac
    [ -z "$best" ] && best="N/A"

    echo "$name,$best,$elapsed,$TIMEOUT,$solved,$note" >> "$OUT_CSV"
    echo "[$name] best=$best time=${elapsed}ms solved=$solved $note"
}

echo "solver,best_cost,time_ms,timeout_s,solved,note" > "$OUT_CSV"
echo ">>> bench_maxsat : $INPUT (timeout ${TIMEOUT}s)"
echo ""

SEED="${SEED:-42}"

run_solver "evalmaxsat" "$EVALMAXSAT_BIN" --TCT "$TIMEOUT" @
run_solver "nuwls"      "$NUWLS_BIN"      @ "$SEED"
run_solver "satlike"    "$SATLIKE_BIN"    @ "$SEED"
run_solver "loandra"    "$LOANDRA_BIN"    @

echo ""
echo ">>> Resultats : $OUT_CSV"
column -t -s ',' "$OUT_CSV"
