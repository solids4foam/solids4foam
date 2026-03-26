#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# cylindricalPressureVessel regression test
# Checks the final point-displacement magnitude at the inner
# radius probe used by the tutorial.
# ============================================================

DISP_MIN=3.10
DISP_MAX=3.25

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "cylindricalPressureVessel regression test"
echo "Final probe displacement magnitude in [${DISP_MIN}, ${DISP_MAX}]"
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
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

value_file=$(find "${CASE_DIR}/postProcessing" -name 'solidPointDisplacement_disp.dat' -print | tail -n 1)
if [[ -z "${value_file}" ]]; then
    echo "FAIL: Could not find point displacement output"
    exit 1
fi

disp_mag=$(awk 'END {print $5}' "${value_file}")
if [[ -z "${disp_mag}" ]]; then
    echo "FAIL: Could not extract probe displacement magnitude"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${disp_mag} >= ${DISP_MIN} && ${disp_mag} <= ${DISP_MAX})}"; then
    printf "PASS: Final probe displacement magnitude = %.6g\n" "${disp_mag}"
else
    printf "FAIL: Final probe displacement magnitude = %.6g\n" "${disp_mag}"
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
