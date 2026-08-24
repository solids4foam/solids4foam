#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# cooksMembrane regression test
# Checks the top-right corner displacement in the elastoplastic case.
# ============================================================

TOP_RIGHT_MIN=0.0068
TOP_RIGHT_MAX=0.0074

ALLRUN_LOGFILE="log.Allrun"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"

echo "============================================================"
echo "cooksMembrane regression test"
echo "Top-right displacement in [${TOP_RIGHT_MIN}, ${TOP_RIGHT_MAX}] m"
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

if [[ ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    echo "FAIL: Could not find ${DISP_FILE}"
    exit 1
fi

top_right_disp=$(awk 'END {print $5}' "${CASE_DIR}/${DISP_FILE}")

if [[ -z "${top_right_disp}" ]]; then
    echo "FAIL: Could not extract top-right displacement"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${top_right_disp} >= ${TOP_RIGHT_MIN} && ${top_right_disp} <= ${TOP_RIGHT_MAX})}"; then
    printf "PASS: Top-right displacement = %.6g\n" "${top_right_disp}"
else
    printf "FAIL: Top-right displacement = %.6g\n" "${top_right_disp}"
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
