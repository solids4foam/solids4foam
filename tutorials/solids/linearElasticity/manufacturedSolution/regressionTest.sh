#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# L2 error tolerances for the coarse hex and tet meshes (OpenFOAM.com v2412).
HEX_DISP_TOL=8e-8
HEX_STRESS_TOL=6e5

TET_DISP_TOL=1.5e-7
TET_STRESS_TOL=9e5

HIGH_ORDER_DISP_TOL=7e-9
HIGH_ORDER_STRESS_TOL=7e4

# Mixed displacement-pressure approach (petscSnesPressure), coarse hex mesh.
# Measured on the 5^3 hex mesh (OpenFOAM.com v2312, PETSc 3.22, macOS):
# displacement 7.94e-8 m, stress 5.26e5 Pa, pressure 2.34e5 Pa (cell RMS).
# The tolerances are about 25% above these values, comparable with the margins
# of the displacement-only tolerances above, for platform and OpenFOAM-variant
# differences; a tighter solve changes these errors by less than 1e-6
# relative. One mesh refinement reduces each error by a factor of 2.7-3.5, so
# a regression of that size is detected
MIXED_HEX_DISP_TOL=1.0e-7
MIXED_HEX_STRESS_TOL=6.6e5
MIXED_HEX_PRESSURE_TOL=2.9e5

APPROACHES=(
    segregated
    petscSnes
    petscSnesPressure
    highOrder-movingLeastSquares
    highOrder-kExactLeastSquares
)

MESHES=(hex tet)

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-}"
: "${FOAM_LD_LIBRARY_PATH:=}"
source solids4FoamScripts.sh

# Usage: ./regressionTest.sh [--check-only] [--mesh hex|tet] [approach ...]
# With approach names, only those approaches are run and checked
CHECK_ONLY=false
SELECTED=()
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
    case "${args[i]}" in
        --check-only|--no-run) CHECK_ONLY=true ;;
        --mesh) i=$((i + 1)); MESHES=("${args[i]}") ;;
        -*) echo "Unknown option ${args[i]}"; exit 1 ;;
        *) SELECTED+=("${args[i]}") ;;
    esac
done
if ((${#SELECTED[@]})); then
    for approach in "${SELECTED[@]}"; do
        known=false
        for candidate in "${APPROACHES[@]}"; do
            [[ "$candidate" == "$approach" ]] && known=true
        done
        if [[ "$known" == false ]]; then
            echo "Unknown approach $approach"
            exit 1
        fi
    done
    APPROACHES=("${SELECTED[@]}")
fi

extract_norm()
{
    local marker="$1"
    local column="$2"
    { grep -A2 "$marker" "$SOLVER_LOG" || true; } | tail -n 1 | awk -v col="$column" '{print $col}'
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
        if [[ "$approach" == "petscSnesPressure" && "$mesh" != "hex" ]]; then
            echo "SKIP: $approach tolerances are calibrated for the hex mesh only"
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
        if [[ "$approach" == "petscSnesPressure" ]]; then
            # Failed linear solves are allowed inside a converged Newton solve
            # (snes_max_linear_solve_fail); the nonlinear outcome decides
            if grep -q 'Nonlinear solve did not converge' "$SOLVER_LOG" \
            || ! grep -q 'Nonlinear solve converged' "$SOLVER_LOG"; then
                echo "FAIL: $approach did not converge"
                failures=$((failures + 1))
                continue
            fi
            # Report failed linear solves by their KSP convergence reason;
            # DIVERGED_ITS is the iteration limit (ksp_max_it)
            n_linear_fail=$(grep -c 'Linear solve did not converge' "$SOLVER_LOG" || true)
            linear_fail_reasons=$(
                grep -o 'Linear solve did not converge due to [A-Z_]*' "$SOLVER_LOG" \
              | awk '{print $NF}' | sort | uniq -c \
              | awk '{printf "%s%s %d", s, $2, $1; s = ", "}' \
              || true
            )
            echo "NOTE: $approach: $n_linear_fail failed linear solve(s)${linear_fail_reasons:+ ($linear_fail_reasons)}"
        elif grep -q 'DIVERGED_' "$SOLVER_LOG"; then
            echo "FAIL: $approach did not converge"
            failures=$((failures + 1))
            continue
        fi

        displacement_l2=$(extract_norm "Writing DDifference field" 3)
        stress_l2=$(extract_norm "Writing sigmaDifference field" 3)

        echo
        echo "Checking $approach ($mesh)"
        if [[ "$approach" == "petscSnesPressure" ]]; then
            # The nonlinear solve must report convergence, and the pressure
            # error must be present: a missing value fails check_tolerance
            pressure_rms=$(extract_norm "Writing pDifference field" 3)
            check_tolerance "displacement L2" "$displacement_l2" "$MIXED_HEX_DISP_TOL" || failures=$((failures + 1))
            check_tolerance "stress L2" "$stress_l2" "$MIXED_HEX_STRESS_TOL" || failures=$((failures + 1))
            check_tolerance "pressure cell RMS" "$pressure_rms" "$MIXED_HEX_PRESSURE_TOL" || failures=$((failures + 1))
        elif [[ "$approach" == highOrder-* ]]; then
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
