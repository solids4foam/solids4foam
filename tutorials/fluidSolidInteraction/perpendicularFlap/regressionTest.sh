#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# perpendicularFlap FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

DISP_MAX_TOL=1e-2      # max displacement absolute tolerance

# Regression end time for the copied case only
REG_END_TIME=4

# Reference values
REF_MAX_DISP=0.118329

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_displacement.dat"

echo "============================================================"
echo "perpendicularFlap FSI regression test"
echo "Max displacement difference < ${DISP_MAX_TOL}"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done

    sed -i "s/^\(endTime[[:space:]]*\).*/\1${REG_END_TIME};/" "${CASE_DIR}/system/controlDict"
}

# ------------------------------------------------------------
# Clean & run case
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

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_displacement() {
    awk '{print $2}' "${CASE_DIR}/${DISP_FILE}" | sort -g | tail -1
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

max_disp=$(extract_max_displacement)

if [[ -z "${max_disp}" ]]; then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

disp_diff=$(awk "BEGIN {print ${max_disp} - ${REF_MAX_DISP}}")
disp_diff_abs=$(abs "${disp_diff}")

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${disp_diff_abs} < ${DISP_MAX_TOL})}"; then
    printf "PASS: max displacement = %.6g (Δ = %.3g)\n" \
        "${max_disp}" "${disp_diff_abs}"
else
    printf "FAIL: max displacement = %.6g (Δ = %.3g)\n" \
        "${max_disp}" "${disp_diff_abs}"
    failures=$((failures + 1))
fi

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
fi

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
