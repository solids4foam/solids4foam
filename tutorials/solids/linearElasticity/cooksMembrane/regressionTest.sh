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
# cooksMembrane regression test
# Checks the vertical displacement at the top-right corner.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

TOP_RIGHT_MIN=0.031
TOP_RIGHT_MAX=0.033

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

APPROACHES=(
    segregated
    petscSnes
    highOrder
)

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

extract_top_right_disp() {
    local value_file
    value_file=$(find "${CASE_DIR}/postProcessing" -name 'solidPointDisplacement_pointDisp.dat' -print | tail -n 1)
    if [[ -z "${value_file}" ]]; then
        echo "FAIL: Could not find point displacement output"
        return 1
    fi

    awk 'END {print $3}' "${value_file}"
}

check_top_right_disp() {
    local approach="$1"
    local top_right_disp

    top_right_disp=$(extract_top_right_disp) || return 1
    if [[ -z "${top_right_disp}" ]]; then
        echo "FAIL: Could not extract top-right displacement for ${approach}"
        return 1
    fi

    if awk "BEGIN {exit !(${top_right_disp} >= ${TOP_RIGHT_MIN} && ${top_right_disp} <= ${TOP_RIGHT_MAX})}"; then
        printf "PASS: Top-right displacement = %.6g\n" "${top_right_disp}"
        return 0
    fi

    printf "FAIL: Top-right displacement = %.6g\n" "${top_right_disp}"
    return 1
}

failures=0

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    for approach in "${APPROACHES[@]}"; do
        echo
        echo "------------------------------------------------------------"
        echo "Testing approach: ${approach}"
        echo "------------------------------------------------------------"

        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
        ( cd "${CASE_DIR}" && ./Allrun "${approach}" > "${ALLRUN_LOGFILE}" 2>&1 )

        if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
            echo "Skipping ${approach} because it is unavailable in this environment"
            continue
        fi

        if ! check_top_right_disp "${approach}"; then
            failures=$((failures + 1))
        fi
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    if ! check_top_right_disp "check-only"; then
        failures=$((failures + 1))
    fi
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
