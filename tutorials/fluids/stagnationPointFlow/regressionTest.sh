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
# stagnationPointFlow regression test
# Checks the velocity and pressure errors against the exact unsteady
# stagnation-point flow at t = 1 on the 20 x 20 mesh, with the outlet in
# fluxConsistentPatches (fixed pressure, fixedGradient velocity) and no
# momentum predictor. The reference norms are the largest and the mean
# magnitude of the errors over the cells. Without the outlet in
# fluxConsistentPatches they are about 2x (velocity) and 5x (pressure, mean)
# larger, so a change that loses the flux-consistent outlet fails
# ============================================================

# ------------------------------------------------------------
# Regression references and tolerance (relative to the largest error)
# ------------------------------------------------------------

REF_UERROR_MAX=0.0408680818996
REF_UERROR_MEAN=0.00384531878240402
REF_PERROR_MAX=0.0253796620659
REF_PERROR_MEAN=0.0134444670649526
REL_TOL=1e-3

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"

echo "============================================================"
echo "stagnationPointFlow regression test"
echo "Velocity and pressure errors at t = 1 match the reference to"
echo "${REL_TOL} of the largest error"
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
# Checks
# ------------------------------------------------------------

if [[ ! -f "${CASE_DIR}/${SOLVER_LOGFILE}" ]]; then
    echo "FAIL: Could not find ${SOLVER_LOGFILE}"
    exit 1
fi

if ! grep -q "Flux-consistent pressure/velocity patches" \
    "${CASE_DIR}/${SOLVER_LOGFILE}"
then
    echo "FAIL: the case did not set up fluxConsistentPatches"
    exit 1
fi

failures=0

latest_time=$(solids4Foam::latestTime "${CASE_DIR}")

if [[ -z "${latest_time}" ]] \
    || ! awk "BEGIN {exit !((${latest_time} - 1)^2 <= 1e-20)}"
then
    echo "FAIL: the case stopped at '${latest_time}', not at t = 1"
    exit 1
fi

solids4Foam::checkFieldNorms "UError at t = 1" \
    "${CASE_DIR}/${latest_time}/UError" \
    "${REF_UERROR_MAX}" "${REF_UERROR_MEAN}" "${REL_TOL}" \
    || failures=$((failures + 1))

solids4Foam::checkFieldNorms "pError at t = 1" \
    "${CASE_DIR}/${latest_time}/pError" \
    "${REF_PERROR_MAX}" "${REF_PERROR_MEAN}" "${REL_TOL}" \
    || failures=$((failures + 1))

# ------------------------------------------------------------
# An outlet velocity condition that neither fixes the value nor is a
# fixedGradient, here inletOutlet, is rejected at start-up
# ------------------------------------------------------------

REJECT_DIR="${REGRESSION_ROOT}/inletOutlet"

if [ "$CHECK_ONLY" = false ]; then
    rm -rf "${REJECT_DIR}"
    cp -a "${CASE_DIR}" "${REJECT_DIR}"
    (
        cd "${REJECT_DIR}"
        ./Allclean > /dev/null 2>&1 || true
        foamDictionary -entry boundaryField/outlet/type -set inletOutlet \
            0/U > /dev/null 2>&1
        foamDictionary -entry boundaryField/outlet/inletValue \
            -add "uniform (0 0 0)" 0/U > /dev/null 2>&1
        blockMesh > log.blockMesh 2>&1
        ulimit -c 0
        solids4Foam > "${SOLVER_LOGFILE}" 2>&1 || true
    ) 2> /dev/null # the FatalError aborts: drop the shell's abort message
fi

if [ "$CHECK_ONLY" = true ]; then
    echo "Check-only mode: skipping the inletOutlet rejection check"
elif grep -q "fluxConsistentPatches: patch outlet needs a fixed pressure" \
    "${REJECT_DIR}/${SOLVER_LOGFILE}" 2>/dev/null
then
    echo "PASS: an inletOutlet outlet velocity is rejected"
else
    echo "FAIL: an inletOutlet outlet velocity was not rejected"
    failures=$((failures + 1))
fi

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    rm -rf "${REJECT_DIR}"
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
