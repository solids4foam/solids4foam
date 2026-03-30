#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# blockInTreacle FSI regression test
# ============================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"
CASE_DIR="${REGRESSION_ROOT}/main"

# ------------------------------------------------------------
# Regression settings
# ------------------------------------------------------------

# Reference values at t=2 (end of simulation)
REF_FINAL_DX=0.1989     # final x-displacement
REF_FINAL_FX=15.77      # final total x-force

# Tolerances (~1%)
DX_TOL=0.002
FX_TOL=0.20

# Log file
ALLRUN_LOGFILE="log.Allrun"

# Data files (relative to case dir)
DISP_FILE="postProcessing/0/solidPointDisplacement_displacement.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

# ------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------

CHECK_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run)
            CHECK_ONLY=true
            ;;
        *)
            ;;
    esac
done

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

abs() { awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'; }

extract_final_dx() {
    local case_dir="$1"
    awk '$1+0==$1 {val=$2} END {print val}' "${case_dir}/${DISP_FILE}"
}

extract_final_fx() {
    local case_dir="$1"
    awk '$1+0==$1 {val=$2} END {print val}' "${case_dir}/${FORCE_FILE}"
}

check_results() {
    local case_dir="$1"
    local failures=0

    local final_dx
    final_dx=$(extract_final_dx "${case_dir}")

    local final_fx
    final_fx=$(extract_final_fx "${case_dir}")

    if [[ -z "${final_dx}" || -z "${final_fx}" ]]; then
        echo "FAIL: Could not extract regression quantities"
        return 1
    fi

    local dx_diff
    dx_diff=$(awk "BEGIN {print ${final_dx} - ${REF_FINAL_DX}}")
    local dx_diff_abs
    dx_diff_abs=$(abs "${dx_diff}")

    local fx_diff
    fx_diff=$(awk "BEGIN {print ${final_fx} - ${REF_FINAL_FX}}")
    local fx_diff_abs
    fx_diff_abs=$(abs "${fx_diff}")

    if awk "BEGIN {exit !(${dx_diff_abs} < ${DX_TOL})}"; then
        printf "PASS: final Dx = %.6g (Δ = %.3g)\n" "${final_dx}" "${dx_diff_abs}"
    else
        printf "FAIL: final Dx = %.6g (Δ = %.3g)\n" "${final_dx}" "${dx_diff_abs}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${fx_diff_abs} < ${FX_TOL})}"; then
        printf "PASS: final Fx = %.6g (Δ = %.3g)\n" "${final_fx}" "${fx_diff_abs}"
    else
        printf "FAIL: final Fx = %.6g (Δ = %.3g)\n" "${final_fx}" "${fx_diff_abs}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

echo "============================================================"
echo "blockInTreacle FSI regression test"
echo "Final Dx difference        < ${DX_TOL}"
echo "Final Fx difference        < ${FX_TOL}"
echo "============================================================"
echo

CASE_DIR="${REGRESSION_ROOT}/main"

if [ "${CHECK_ONLY}" = false ]; then
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        [[ "${base_item}" == "regressionTests" ]] && continue
        cp -a "${item}" "${CASE_DIR}/"
    done

    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

    # Run case
    ( cd "${CASE_DIR}" && ./Allrun ) > "${CASE_DIR}/${ALLRUN_LOGFILE}" 2>&1

    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
    CASE_DIR="${SCRIPT_DIR}"
fi

failures=0
check_results "${CASE_DIR}" || failures=$((failures + $?))

echo
if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED"
    echo "============================================================"
    exit 0
else
    echo "============================================================"
    echo "Regression test FAILED (${failures} checks)"
    echo "============================================================"
    exit 1
fi
