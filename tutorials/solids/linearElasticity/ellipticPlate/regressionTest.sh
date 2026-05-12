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
# ellipticPlate regression test
# Checks selected solution approaches against the expected
# solution extrema on the default mesh.
# ============================================================

EPS_MIN=3.8e-5
EPS_MAX=4.1e-5
SIGMA_MIN=9.65e6
SIGMA_MAX=9.75e6

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

APPROACHES=(
    segregated
    petscSnes
)

echo "============================================================"
echo "ellipticPlate regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
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
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null \
        | tail -n 1 \
        | awk '{print $NF}' || true
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null \
        | tail -n 1 \
        | awk '{print $NF}' || true
}

check_solver_extrema() {
    local approach="$1"
    local epsilon
    local sigma
    local failures=0

    epsilon=$(extract_max_epsilon)
    sigma=$(extract_max_sigma)

    if [[ -z "${epsilon}" || -z "${sigma}" ]]; then
        echo "FAIL: Could not extract one or more regression quantities for ${approach}"
        return 1
    fi

    if awk "BEGIN {exit !(${epsilon} >= ${EPS_MIN} && ${epsilon} <= ${EPS_MAX})}"; then
        printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
    else
        printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${sigma} <= ${SIGMA_MAX})}"; then
        printf "PASS: Max sigmaEq = %.6g\n" "${sigma}"
    else
        printf "FAIL: Max sigmaEq = %.6g\n" "${sigma}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

failures=0

if [ "$CHECK_ONLY" = false ]; then
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

        if ! check_solver_extrema "${approach}"; then
            failures=$((failures + 1))
        fi
    done
else
    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    if ! check_solver_extrema "check-only"; then
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
