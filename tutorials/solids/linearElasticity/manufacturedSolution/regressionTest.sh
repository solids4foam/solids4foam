#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# L2 error tolerances for the coarse hex and tet meshes (OpenFOAM.com v2412).
HEX_DISP_TOL=8e-8
HEX_STRESS_TOL=6e5

TET_DISP_TOL=1.5e-7
TET_STRESS_TOL=9e5

HIGH_ORDER_DISP_TOL=5e-9
HIGH_ORDER_STRESS_TOL=7e4

APPROACHES=(
    segregated
    petscSnes
    highOrder-movingLeastSquares
    highOrder-kExactLeastSquares
)

MESHES=(hex tet)

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

check_log_message()
{
    local message="$1"

    if ! grep -Fq "$message" "$SOLVER_LOG"; then
        echo "FAIL: expected '$message' in $SOLVER_LOG"
        return 1
    fi
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

    if awk "BEGIN {exit !($value <= $tolerance)}"; then
        printf 'PASS: %s = %.8g\n' "$label" "$value"
    else
        printf 'FAIL: %s = %.8g, expected <= %g\n' \
            "$label" "$value" "$tolerance"
        return 1
    fi
}

failures=0
for mesh in "${MESHES[@]}"; do
    if [[ "$mesh" == "tet" ]] && ! command -v gmsh >/dev/null 2>&1; then
        echo "SKIP: tet mesh requires Gmsh"
        continue
    fi

    for approach in "${APPROACHES[@]}"; do
        if [[ "$approach" != "segregated" && -z "${PETSC_DIR:-}" ]]; then
            echo "SKIP: $approach requires PETSc"
            continue
        fi

        CASE_DIR="$SCRIPT_DIR/regressionTests/$approach-$mesh"

        if [[ "$CHECK_ONLY" == false ]]; then
            rm -rf "$CASE_DIR"
            mkdir -p "$CASE_DIR"
            for item in "$SCRIPT_DIR"/*; do
                name=$(basename "$item")
                case "$name" in
                    regressionTests|verification|postProcessing|processor*|log.*|[1-9]*) continue ;;
                esac
                cp -a "$item" "$CASE_DIR/"
            done
            rm -rf "$CASE_DIR/constant/polyMesh"
            (cd "$CASE_DIR" && ./Allrun "$approach" "$mesh" > log.Allrun 2>&1)
        fi

        SOLVER_LOG="$CASE_DIR/log.solids4Foam"
        if [[ ! -f "$SOLVER_LOG" ]] || ! grep -q '^End' "$SOLVER_LOG"; then
            echo "FAIL: missing or incomplete solver log: $SOLVER_LOG"
            failures=$((failures + 1))
            continue
        fi
        if grep -q 'DIVERGED_' "$SOLVER_LOG"; then
            echo "FAIL: $approach did not converge"
            failures=$((failures + 1))
            continue
        fi

        displacement_l2=$(extract_norm "Writing DDifference field" 3)
        stress_l2=$(extract_norm "Writing sigmaDifference field" 3)

        echo
        echo "Checking $approach ($mesh)"
        if [[ "$approach" == highOrder-* ]]; then
            check_log_message 'Using volume-averaged manufactured body force' || failures=$((failures + 1))
            if [[ "$approach" == "highOrder-kExactLeastSquares" ]]; then
                check_log_message 'Using cell-average analytical displacement' || failures=$((failures + 1))
            else
                check_log_message 'Using point-valued analytical displacement' || failures=$((failures + 1))
            fi
            check_tolerance "displacement L2" "$displacement_l2" "$HIGH_ORDER_DISP_TOL" || failures=$((failures + 1))
            check_tolerance "stress L2" "$stress_l2" "$HIGH_ORDER_STRESS_TOL" || failures=$((failures + 1))
        elif [[ "$mesh" == "tet" ]]; then
            check_tolerance "displacement L2" "$displacement_l2" "$TET_DISP_TOL" || failures=$((failures + 1))
            check_tolerance "stress L2" "$stress_l2" "$TET_STRESS_TOL" || failures=$((failures + 1))
        else
            check_tolerance "displacement L2" "$displacement_l2" "$HEX_DISP_TOL" || failures=$((failures + 1))
            check_tolerance "stress L2" "$stress_l2" "$HEX_STRESS_TOL" || failures=$((failures + 1))
        fi
    done
done

if ((failures)); then
    echo "Regression test FAILED ($failures checks)"
    exit 1
fi

echo "Regression test PASSED"
