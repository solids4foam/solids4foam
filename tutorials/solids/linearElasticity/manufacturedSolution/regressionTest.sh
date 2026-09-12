#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CASE_DIR="${SCRIPT_DIR}/regressionTests/main"

export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-}"
: "${FOAM_LD_LIBRARY_PATH:=}"
source solids4FoamScripts.sh
solids4Foam::caseDoesNotRunWithFoamExtend
solids4Foam::caseDoesNotRunWithOpenFOAMOrg

CHECK_ONLY=false
for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run) CHECK_ONLY=true ;;
    esac
done

if [[ "$CHECK_ONLY" == false ]]; then
    rm -rf "$CASE_DIR"
    mkdir -p "$CASE_DIR"
    for item in "$SCRIPT_DIR"/*; do
        name=$(basename "$item")
        case "$name" in
            regressionTests|verification|log.*|[1-9]*) continue ;;
        esac
        cp -a "$item" "$CASE_DIR/"
    done
    rm -rf "$CASE_DIR/constant/polyMesh"
    (cd "$CASE_DIR" && ./Allrun > log.Allrun 2>&1)
fi

SOLVER_LOG="$CASE_DIR/log.solids4Foam"
if [[ ! -f "$SOLVER_LOG" ]]; then
    echo "FAIL: solver log not found: $SOLVER_LOG"
    exit 1
fi

extract_norm()
{
    local marker="$1"
    local column="$2"
    grep -A2 "$marker" "$SOLVER_LOG" | tail -n 1 | awk -v col="$column" '{print $col}'
}

displacement_l2=$(extract_norm "Writing DDifference field" 3)
displacement_linf=$(extract_norm "Writing DDifference field" 4)
stress_l2=$(extract_norm "Writing sigmaDifference field" 3)
stress_linf=$(extract_norm "Writing sigmaDifference field" 4)

check_range()
{
    local label="$1"
    local value="$2"
    local minimum="$3"
    local maximum="$4"

    if awk "BEGIN {exit !($value >= $minimum && $value <= $maximum)}"; then
        printf 'PASS: %s = %.8g\n' "$label" "$value"
    else
        printf 'FAIL: %s = %.8g, expected [%g, %g]\n' \
            "$label" "$value" "$minimum" "$maximum"
        return 1
    fi
}

failures=0
check_range "displacement L2" "$displacement_l2" 4e-8 8e-8 || failures=$((failures + 1))
check_range "displacement Linf" "$displacement_linf" 1e-7 1.5e-7 || failures=$((failures + 1))
check_range "stress L2" "$stress_l2" 4e5 6e5 || failures=$((failures + 1))
check_range "stress Linf" "$stress_linf" 1.4e6 1.9e6 || failures=$((failures + 1))

if [[ "$CHECK_ONLY" == false ]]; then
    (cd "$CASE_DIR" && ./Allclean >/dev/null 2>&1) || true
fi

if ((failures)); then
    echo "Regression test FAILED ($failures checks)"
    exit 1
fi

echo "Regression test PASSED"
