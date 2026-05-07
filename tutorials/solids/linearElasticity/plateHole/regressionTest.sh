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
# Plate-with-hole regression test
# Checks selected solution approaches against the analytical solution
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

DISP_TOL=1e-7
POINT_DISP_TOL=1e-7
STRESS_TOL=2e5   # component-0 LInf

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

APPROACHES=(
    segregated
    petscSnes
    highOrder
)

echo "============================================================"
echo "Plate-with-hole regression test"
echo "DDifference LInf        < ${DISP_TOL}"
echo "pointDDifference LInf   < ${POINT_DISP_TOL}"
echo "Stress component-0 LInf < ${STRESS_TOL}"
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
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_disp_linf() {
    local field="$1"
    grep -A2 "Writing ${field} field" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | grep "Norms:" -A1 \
        | tail -n 1 \
        | awk '{print $3}' || true
}

extract_stress_linf_comp0() {
    grep -A6 "Writing cellStressDifference field" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -6 \
        | awk '
            /Component:[[:space:]]*0/ {getline; getline; print $3}
        ' || true
}

check_analytical_norms() {
    local approach="$1"
    local disp_linf
    local point_disp_linf
    local stress_linf
    local failures=0

    disp_linf=$(extract_disp_linf "DDifference")
    point_disp_linf=$(extract_disp_linf "pointDDifference")
    stress_linf=$(extract_stress_linf_comp0)

    if [[ -z "${disp_linf}" || -z "${point_disp_linf}" || -z "${stress_linf}" ]]; then
        echo "FAIL: Could not extract one or more error norms for ${approach}"
        return 1
    fi

    if awk "BEGIN {exit !(${disp_linf} < ${DISP_TOL})}"; then
        printf "PASS: DDifference LInf = %.6g\n" "${disp_linf}"
    else
        printf "FAIL: DDifference LInf = %.6g\n" "${disp_linf}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${point_disp_linf} < ${POINT_DISP_TOL})}"; then
        printf "PASS: pointDDifference LInf = %.6g\n" "${point_disp_linf}"
    else
        printf "FAIL: pointDDifference LInf = %.6g\n" "${point_disp_linf}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${stress_linf} < ${STRESS_TOL})}"; then
        printf "PASS: stress component-0 LInf = %.6g\n" "${stress_linf}"
    else
        printf "FAIL: stress component-0 LInf = %.6g\n" "${stress_linf}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

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

        if ! check_analytical_norms "${approach}"; then
            failures=$((failures + 1))
        fi
    done
else
    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    if ! check_analytical_norms "check-only"; then
        failures=$((failures + 1))
    fi
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
