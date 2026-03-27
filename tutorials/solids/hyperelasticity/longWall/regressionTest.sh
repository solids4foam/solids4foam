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
# longWall regression test
# Uses top-surface stress and displacement histories.
# ============================================================

UY_MIN=0.405
UY_MAX=0.407
SYY_MIN=9.99e7
SYY_MAX=1.001e8

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "longWall regression test"
echo "Top-surface uy in [${UY_MIN}, ${UY_MAX}] m"
echo "Top-surface sigma_yy in [${SYY_MIN}, ${SYY_MAX}] Pa"
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

find_history_file() {
    local name="$1"
    find "${CASE_DIR}/postProcessing" -name "${name}" -print 2>/dev/null | tail -n 1
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

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

disp_file=$(find_history_file 'solidDisplacementstop.dat')
stress_file=$(find_history_file 'solidStressestop.dat')

if [[ -z "${disp_file}" || -z "${stress_file}" ]]; then
    echo "FAIL: Could not find one or more history files"
    exit 1
fi

top_uy=$(awk 'END {print $9}' "${disp_file}")
top_syy=$(awk 'END {print $5}' "${stress_file}")

if [[ -z "${top_uy}" || -z "${top_syy}" ]]; then
    echo "FAIL: Could not extract top-surface values"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${top_uy} >= ${UY_MIN} && ${top_uy} <= ${UY_MAX})}"; then
    printf "PASS: Top-surface uy = %.6g\n" "${top_uy}"
else
    printf "FAIL: Top-surface uy = %.6g\n" "${top_uy}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${top_syy} >= ${SYY_MIN} && ${top_syy} <= ${SYY_MAX})}"; then
    printf "PASS: Top-surface sigma_yy = %.6g\n" "${top_syy}"
else
    printf "FAIL: Top-surface sigma_yy = %.6g\n" "${top_syy}"
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
