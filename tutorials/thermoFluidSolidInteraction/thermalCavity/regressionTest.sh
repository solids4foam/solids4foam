#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# thermalCavity regression test
# Uses the fsiConvergenceData output as a cheap convergence check.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

N_FSI_CORRECTORS_MAX=50

# The final nFsiCorrectors of the removed legacy mechanicalModel, from the last
# commit that had it (mcl-stage8-coverage, c3a92b3d), the same on every fork.
# thermalSolid asks the constitutive implementation only for the density, so
# the same rho drives the same conjugate heat transfer and the same coupling
# iteration count: the framework reproduced it exactly, and must
LEGACY_N_FSI_CORRECTORS=3

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "thermalCavity regression test"
echo "nFsiCorrectors < ${N_FSI_CORRECTORS_MAX}"
echo "============================================================"
echo

prepare_case() {
    local d="${1:-${CASE_DIR}}"
    rm -rf "${d}"
    mkdir -p "${d}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${d}/"
    done
}

shorten_case() {
    local controlDict="${1:-${CASE_DIR}}/system/controlDict"
    sed -i.bak 's/^endTime[[:space:]]\+10;/endTime         0.1;/' "${controlDict}"
    rm -f "${controlDict}.bak"
}

find_fsi_data() {
    local root="${1:-${CASE_DIR}}"
    local candidate
    for candidate in \
        "${root}/postProcessing/0/fsiConvergenceData.dat" \
        "${root}/postProcessing/fluid/0/fsiConvergenceData.dat" \
        "${root}/postProcessing/solid/0/fsiConvergenceData.dat"
    do
        if [[ -f "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done

    find "${root}/postProcessing" -name 'fsiConvergenceData.dat' \
        -print 2>/dev/null | tail -n 1
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
    prepare_case "${CASE_DIR}"
    shorten_case "${CASE_DIR}"
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

fsi_data=$(find_fsi_data)
if [[ -z "${fsi_data}" ]]; then
    echo "FAIL: Could not find fsiConvergenceData output"
    exit 1
fi

n_fsi_correctors=$(grep -v '^[[:space:]]*#' "${fsi_data}" | tail -n 1 | awk '{print $2}')

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

# ------------------------------------------------------------
# Against the legacy answer
# ------------------------------------------------------------

main_log=$(find "${CASE_DIR}" -name "${SOLVER_LOGFILE}" | tail -n 1)

if [[ -n "${main_log}" ]] \
    && grep -q "Selecting mechanical constitutive law" "${main_log}"
then
    echo "PASS: the solid took its density from the framework"
else
    echo "FAIL: the solver log shows no mechanical constitutive law"
    failures=$((failures + 1))
fi

if [[ "${n_fsi_correctors}" == "${LEGACY_N_FSI_CORRECTORS}" ]]; then
    printf "PASS: nFsiCorrectors = %s, as the legacy model\n" \
        "${n_fsi_correctors}"
else
    printf "FAIL: nFsiCorrectors differ (%s, legacy model %s)\n" \
        "${n_fsi_correctors}" "${LEGACY_N_FSI_CORRECTORS}"
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
