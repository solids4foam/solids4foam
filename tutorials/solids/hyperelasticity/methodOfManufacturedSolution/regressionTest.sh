#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Regression test for the methodOfManufacturedSolution tutorial: runs the
# default coarse steady and transient cases with the segregated and PETSc SNES
# approaches and checks the final displacement and stress L2 errors against
# the manufactured solution (tolerances from OpenFOAM.com v2412).

STEADY_DISP_TOL=1.5e-3
STEADY_STRESS_TOL=1.8e4

TRANSIENT_DISP_TOL=2.0e-3
TRANSIENT_STRESS_TOL=1.9e4

APPROACHES=(segregated petscSnes)
CASES=(steady transient)

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
for case_name in "${CASES[@]}"; do
    for approach in "${APPROACHES[@]}"; do
        if [[ "$approach" == "petscSnes" && -z "${PETSC_DIR:-}" ]]; then
            echo "SKIP: $approach requires PETSc"
            continue
        fi

        RUN_DIR="$SCRIPT_DIR/regressionTests/$case_name-$approach"
        CASE_DIR="$RUN_DIR/$case_name"

        if [[ "$CHECK_ONLY" == false ]]; then
            rm -rf "$RUN_DIR"
            mkdir -p "$RUN_DIR"
            cp -a "$SCRIPT_DIR/src" "$SCRIPT_DIR/$case_name" "$RUN_DIR/"
            rm -rf "$RUN_DIR/src/lnInclude" "$CASE_DIR/constant/polyMesh"
            (cd "$CASE_DIR" && ./Allclean > /dev/null 2>&1 \
                && ./Allrun "$approach" > log.Allrun 2>&1)
        fi

        SOLVER_LOG="$CASE_DIR/log.solids4Foam"
        if [[ ! -f "$SOLVER_LOG" ]] || ! grep -q '^End' "$SOLVER_LOG"; then
            echo "FAIL: missing or incomplete solver log: $SOLVER_LOG"
            failures=$((failures + 1))
            continue
        fi
        if grep -q 'DIVERGED_' "$SOLVER_LOG"; then
            echo "FAIL: $case_name $approach did not converge"
            failures=$((failures + 1))
            continue
        fi

        displacement_l2=$(extract_norm "Writing DDifference field" 3)
        stress_l2=$(extract_norm "Writing sigmaDifference field" 3)

        echo
        echo "Checking $case_name ($approach)"
        if [[ "$case_name" == "steady" ]]; then
            check_tolerance "displacement L2" "$displacement_l2" "$STEADY_DISP_TOL" || failures=$((failures + 1))
            check_tolerance "stress L2" "$stress_l2" "$STEADY_STRESS_TOL" || failures=$((failures + 1))
        else
            check_tolerance "displacement L2" "$displacement_l2" "$TRANSIENT_DISP_TOL" || failures=$((failures + 1))
            check_tolerance "stress L2" "$stress_l2" "$TRANSIENT_STRESS_TOL" || failures=$((failures + 1))
        fi
    done
done

if ((failures)); then
    echo "Regression test FAILED ($failures checks)"
    exit 1
fi

echo "Regression test PASSED"
