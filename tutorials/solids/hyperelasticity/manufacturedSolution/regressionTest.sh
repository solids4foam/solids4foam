#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Regression test for the manufacturedSolution tutorial: runs the default
# coarse case in the steady and dynamic modes with the segregated and PETSc
# SNES approaches and checks the final displacement and stress L2 errors
# against the manufactured solution (tolerances from OpenFOAM.com v2412).

STEADY_DISP_TOL=1.5e-3
STEADY_STRESS_TOL=1.8e4

DYNAMIC_DISP_TOL=2.0e-3
DYNAMIC_STRESS_TOL=1.9e4

# Both high-order reconstructions (cubic) share the same tolerances
STEADY_HIGH_ORDER_DISP_TOL=8e-5
STEADY_HIGH_ORDER_STRESS_TOL=1e3

DYNAMIC_HIGH_ORDER_DISP_TOL=8e-4
DYNAMIC_HIGH_ORDER_STRESS_TOL=4e3

APPROACHES=(
    segregated
    petscSnes
    highOrder-movingLeastSquares
    highOrder-kExactLeastSquares
)
MODES=(steady dynamic)

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-}"
: "${FOAM_LD_LIBRARY_PATH:=}"
source solids4FoamScripts.sh

CHECK_ONLY=false
for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run) CHECK_ONLY=true ;;
    esac
done

check_log_message()
{
    local message="$1"

    if ! grep -Fq "$message" "$SOLVER_LOG"; then
        echo "FAIL: expected '$message' in $SOLVER_LOG"
        return 1
    fi
}

extract_norm()
{
    local marker="$1"
    local column="$2"
    grep -A2 "$marker" "$SOLVER_LOG" | tail -n 1 | awk -v col="$column" '{print $col}'
}

check_tolerance()
{
    local label="$1"
    local value="$2"
    local tolerance="$3"

    if [[ ! "$value" =~ ^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$ ]]; then
        echo "FAIL: $label has missing or invalid value '$value'"
        return 1
    fi

    # Errors must be positive: zero indicates a degenerate comparison
    if awk "BEGIN {exit !($value > 0 && $value <= $tolerance)}"; then
        printf 'PASS: %s = %.8g\n' "$label" "$value"
    else
        printf 'FAIL: %s = %.8g, expected in (0, %g]\n' \
            "$label" "$value" "$tolerance"
        return 1
    fi
}

failures=0
for mode in "${MODES[@]}"; do
    for approach in "${APPROACHES[@]}"; do
        if [[ "$approach" != "segregated" && -z "${PETSC_DIR:-}" ]]; then
            echo "SKIP: $approach requires PETSc"
            continue
        fi

        CASE_DIR="$SCRIPT_DIR/regressionTests/$mode-$approach"

        if [[ "$CHECK_ONLY" == false ]]; then
            rm -rf "$CASE_DIR"
            mkdir -p "$CASE_DIR"
            for item in "$SCRIPT_DIR"/*; do
                name=$(basename "$item")
                case "$name" in
                    regressionTests|verification|docs|postProcessing|processor*|log.*|[1-9]*) continue ;;
                esac
                cp -a "$item" "$CASE_DIR/"
            done
            rm -rf "$CASE_DIR/src/lnInclude" "$CASE_DIR/constant/polyMesh"
            (cd "$CASE_DIR" && ./Allclean > /dev/null 2>&1 \
                && ./Allrun "$mode" "$approach" > log.Allrun 2>&1)
        fi

        SOLVER_LOG="$CASE_DIR/log.solids4Foam"
        if [[ ! -f "$SOLVER_LOG" ]] || ! grep -q '^End' "$SOLVER_LOG"; then
            echo "FAIL: missing or incomplete solver log: $SOLVER_LOG"
            failures=$((failures + 1))
            continue
        fi
        if grep -q 'DIVERGED_' "$SOLVER_LOG"; then
            echo "FAIL: $mode $approach did not converge"
            failures=$((failures + 1))
            continue
        fi

        displacement_l2=$(extract_norm "Writing DDifference field" 3)
        stress_l2=$(extract_norm "Writing sigmaDifference field" 3)

        echo
        echo "Checking $mode ($approach)"
        if [[ "$approach" == highOrder-* ]]; then
            check_log_message 'Using volume-averaged manufactured body force' || failures=$((failures + 1))
            if [[ "$approach" == "highOrder-kExactLeastSquares" ]]; then
                check_log_message 'Using cell-average analytical displacement' || failures=$((failures + 1))
            fi
            if [[ "$mode" == "steady" ]]; then
                disp_tol="$STEADY_HIGH_ORDER_DISP_TOL"
                stress_tol="$STEADY_HIGH_ORDER_STRESS_TOL"
            else
                disp_tol="$DYNAMIC_HIGH_ORDER_DISP_TOL"
                stress_tol="$DYNAMIC_HIGH_ORDER_STRESS_TOL"
            fi
        elif [[ "$mode" == "steady" ]]; then
            disp_tol="$STEADY_DISP_TOL"
            stress_tol="$STEADY_STRESS_TOL"
        else
            disp_tol="$DYNAMIC_DISP_TOL"
            stress_tol="$DYNAMIC_STRESS_TOL"
        fi
        check_tolerance "displacement L2" "$displacement_l2" "$disp_tol" || failures=$((failures + 1))
        check_tolerance "stress L2" "$stress_l2" "$stress_tol" || failures=$((failures + 1))
    done
done

if ((failures)); then
    echo "Regression test FAILED ($failures checks)"
    exit 1
fi

echo "Regression test PASSED"
