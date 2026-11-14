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
# thinPlateInDuct regression test (experimental thinBodyBlocking)
# Runs the 0.43-cell plate with slip walls to t = 8 s and checks the
# time-average (2-8 s) of Q/(w*A) of the outlet flux: leaking with the
# switch off (default), exact with crossedFaces
# ============================================================

OFF_MIN=0.5
OFF_MAX=0.95
CROSSED_MIN=0.999
CROSSED_MAX=1.001

ALLRUN_LOGFILE="log.Allrun"

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"
    for item in "${SCRIPT_DIR}"/*; do
        if [[ "$(basename "${item}")" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done
}

CHECK_ONLY=false
for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run) CHECK_ONLY=true ;;
        *) ;;
    esac
done

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

failures=0
check_case() {
    local name=$1 min=$2 max=$3
    if ( cd "${CASE_DIR}" && python3 checkFlux.py "${name}" "${min}" "${max}" ); then
        echo "PASS: ${name}"
    else
        echo "FAIL: ${name}"
        failures=$((failures + 1))
    fi
}

check_case run_slip_0.43_off "${OFF_MIN}" "${OFF_MAX}"
check_case run_slip_0.43_crossedFaces "${CROSSED_MIN}" "${CROSSED_MAX}"

if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
fi

echo
if (( failures == 0 )); then
    echo "Regression test PASSED"
    exit 0
else
    echo "Regression test FAILED (${failures} checks)"
    exit 1
fi
