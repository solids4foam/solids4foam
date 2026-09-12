#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# cantileverVibration regression test
#
# Checks the peak magnitude of the tip displacement over the run (one full
# oscillation, 0 <= t <= 0.65 s), taken from column 5 (magD) of
# postProcessing/0/solidPointDisplacement_pointDisp.dat.
#
# The peak is used rather than the final-time value because, at the end of
# the period, the tip is close to its undeformed position: the final value is
# small and very sensitive to small phase errors, whereas the peak amplitude is
# a well-conditioned measure of the dynamic, geometrically nonlinear response.
#
# Measured with the default petscSnes approach (6 x 6 x 60 mesh,
# deltaT = 0.005 s, OpenFOAM-v2512): peak = 2.7247 m at t = 0.31 s.
# The Abaqus (C3D8) reference peak is 2.8007 m at t = 0.318 s (see
# abaqusC3D8.dat); the tutorial mesh is a coarse demonstration mesh, so the
# band below is centred on the measured solids4foam value (approx. +/- 0.9%),
# not on the Abaqus value.
# ============================================================

PEAK_MIN=2.70
PEAK_MAX=2.75

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "cantileverVibration regression test"
echo "Peak tip displacement magnitude in [${PEAK_MIN}, ${PEAK_MAX}] m"
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
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

failures=0

# The Allrun script can exit successfully even if the solver fails, so check
# the solver log directly
solver_log="${CASE_DIR}/${SOLVER_LOGFILE}"
if [[ ! -f "${solver_log}" ]]; then
    echo "FAIL: Could not find ${SOLVER_LOGFILE}"
    failures=$((failures + 1))
elif grep -Eq 'FOAM FATAL|^ERROR$|\[stack trace\]' "${solver_log}" \
    || ! grep -q '^End' "${solver_log}"
then
    echo "FAIL: solids4Foam did not complete successfully"
    failures=$((failures + 1))
else
    echo "PASS: solids4Foam completed"
fi

disp_file=$(find "${CASE_DIR}/postProcessing" \
    -name 'solidPointDisplacement_pointDisp.dat' -print 2>/dev/null | tail -n 1)

if [[ -z "${disp_file}" ]]; then
    echo "FAIL: Could not find point displacement output"
    exit 1
fi

peak=$(awk '!/^#/ && NF >= 5 {if (!n++ || $5 > m) {m = $5; t = $1}} END {if (n) print m}' "${disp_file}")
peak_time=$(awk '!/^#/ && NF >= 5 {if (!n++ || $5 > m) {m = $5; t = $1}} END {if (n) print t}' "${disp_file}")

if [[ -z "${peak}" ]]; then
    echo "FAIL: Could not extract the peak tip displacement"
    exit 1
fi

if awk "BEGIN {exit !(${peak} >= ${PEAK_MIN} && ${peak} <= ${PEAK_MAX})}"; then
    printf "PASS: Peak tip displacement = %.6g m at t = %s s\n" "${peak}" "${peak_time}"
else
    printf "FAIL: Peak tip displacement = %.6g m at t = %s s\n" "${peak}" "${peak_time}"
    failures=$((failures + 1))
fi

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
