#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"
SOLIDS4FOAM_ROOT_ABS=$(cd "${SCRIPT_DIR}/../../../../" && pwd)

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# cantilever2d regression test
# Checks selected solution approaches against the analytical
# benchmark output.
# ============================================================

EPS_MIN=3.5e-4
EPS_MAX=5.0e-4
SIGMA_MIN=8.5e7
SIGMA_MAX=1.05e8
HIGH_ORDER_DISP_TOL=1e-11

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

APPROACHES=(
    petscSnes
    highOrder
)

echo "============================================================"
echo "cantilever2d regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "High-order DDifference LInf < ${HIGH_ORDER_DISP_TOL}"
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

    # The regression copy lives deeper than the source tutorial, so the
    # relative SOLIDS4FOAM_ROOT in this local library build no longer points to
    # the repository root.
    sed -i.bak \
        "s|^SOLIDS4FOAM_ROOT := .*|SOLIDS4FOAM_ROOT := ${SOLIDS4FOAM_ROOT_ABS}|" \
        "${CASE_DIR}/src/Make/options"

    if [[ -f "${SCRIPT_DIR}/constant/polyMesh/blockMeshDict" ]] \
        && [[ ! -f "${CASE_DIR}/constant/polyMesh/blockMeshDict" ]]; then
        mkdir -p "${CASE_DIR}/constant/polyMesh"
        cp -a "${SCRIPT_DIR}/constant/polyMesh/blockMeshDict" \
            "${CASE_DIR}/constant/polyMesh/blockMeshDict"
    fi
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

extract_disp_linf() {
    grep -A2 "Writing DDifference field" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | grep "Norms:" -A1 \
        | tail -n 1 \
        | awk '{print $3}' \
        || true
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

check_high_order_errors() {
    local displacement
    local failures=0

    displacement=$(extract_disp_linf)

    if [[ -z "${displacement}" ]]; then
        echo "FAIL: Could not extract high-order DDifference LInf"
        return 1
    fi

    if awk "BEGIN {exit !(${displacement} < ${HIGH_ORDER_DISP_TOL})}"; then
        printf "PASS: High-order DDifference LInf = %.6g\n" "${displacement}"
    else
        printf "FAIL: High-order DDifference LInf = %.6g\n" "${displacement}"
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

        if [[ "${approach}" == "highOrder" ]] && ! check_high_order_errors; then
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

    if grep -q "highOrderResidual true" "${CASE_DIR}/constant/solidProperties" \
        && ! check_high_order_errors; then
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
