#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# thermalCavity regression test
# Uses the fsiConvergenceData output as a cheap convergence check.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

N_FSI_CORRECTORS_MAX=50

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
FSI_DATA="postProcessing/0/fsiConvergenceData.dat"

echo "============================================================"
echo "thermalCavity regression test"
echo "nFsiCorrectors < ${N_FSI_CORRECTORS_MAX}"
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
}

shorten_case() {
    local controlDict="${CASE_DIR}/system/controlDict"
    sed -i.bak 's/^endTime[[:space:]]\+10;/endTime         0.1;/' "${controlDict}"
    rm -f "${controlDict}.bak"
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
    shorten_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    shorten_case
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

if [[ ! -f "${CASE_DIR}/${FSI_DATA}" ]]; then
    echo "FAIL: Could not find ${FSI_DATA}"
    exit 1
fi

n_fsi_correctors=$(grep -v '^[[:space:]]*#' "${CASE_DIR}/${FSI_DATA}" | tail -n 1 | awk '{print $2}')

if [[ -z "${n_fsi_correctors}" ]]; then
    echo "FAIL: Could not extract nFsiCorrectors"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${n_fsi_correctors} < ${N_FSI_CORRECTORS_MAX})}"; then
    printf "PASS: nFsiCorrectors = %.6g\n" "${n_fsi_correctors}"
else
    printf "FAIL: nFsiCorrectors = %.6g\n" "${n_fsi_correctors}"
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
