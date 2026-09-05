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
    highOrder-movingLeastSquares
    highOrder-kExactLeastSquares
)

# Settings for Test-leastSquaresS4fGrad test
GRADIENT_REL_TOL=5e-1
GRADIENT_ABS_TOL=1e-10
GRADIENT_LOGFILE="log.Test-leastSquaresS4fGrad.serial"

echo "============================================================"
echo "cooksMembrane regression test"
echo "Top-right displacement in [${TOP_RIGHT_MIN}, ${TOP_RIGHT_MAX}] m"
echo "Gradient relative tolerance: ${GRADIENT_REL_TOL}"
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

select_run_approach() {
    local requested="$1"

    case "${requested}" in
        highOrder-movingLeastSquares|highOrder-kExactLeastSquares)
            local least_squares_type="${requested#highOrder-}"
            sed -E -i.bak \
                "s/^([[:space:]]*)type[[:space:]]+(movingLeastSquares|kExactLeastSquares);/\\1type ${least_squares_type};/" \
                "${CASE_DIR}/constant/solidProperties.highOrder"
            rm -f "${CASE_DIR}/constant/solidProperties.highOrder.bak"
            RUN_APPROACH=highOrder
            ;;
        *)
            RUN_APPROACH="${requested}"
            ;;
    esac
}

check_gradient_log() {
    local failures=0

    check_gradient_reference_field phi_x 0 0 0 0 \
        || failures=$((failures + 1))
    check_gradient_reference_field phi_linear 0 0 0 0 \
        || failures=$((failures + 1))
    check_gradient_reference_field \
        phi_quadratic 0.000721848 0.00377471 0.00165316 0.00377471 \
        || failures=$((failures + 1))
    check_gradient_reference_field \
        phi_quadratic_cross 0.00034808 0.00216383 0.000796681 0.00216383 \
        || failures=$((failures + 1))
    check_gradient_reference_field \
        phi_trigonometric 0.0150774 0.095765 0.0345717 0.095765 \
        || failures=$((failures + 1))

    return "${failures}"
}

extract_gradient_row() {
    local field="$1"

    awk -v field="${field}" \
        '$1 == field && NF >= 5 {print $2, $3, $4, $5; exit}' \
        "${CASE_DIR}/${GRADIENT_LOGFILE}" 2>/dev/null || true
}

check_gradient_reference_field() {
    local field="$1"
    local ref_l2="$2"
    local ref_linf="$3"
    local ref_boundary_l2="$4"
    local ref_boundary_linf="$5"
    local row
    local l2 linf boundary_l2 boundary_linf

    row=$(extract_gradient_row "${field}")
    if [[ -z "${row}" ]]; then
        echo "FAIL: Could not extract gradient norms for ${field}"
        return 1
    fi

    IFS=' ' read -r l2 linf boundary_l2 boundary_linf <<< "${row}"
    if awk -v l2="${l2}" -v linf="${linf}" \
        -v boundary_l2="${boundary_l2}" \
        -v boundary_linf="${boundary_linf}" \
        -v ref_l2="${ref_l2}" -v ref_linf="${ref_linf}" \
        -v ref_boundary_l2="${ref_boundary_l2}" \
        -v ref_boundary_linf="${ref_boundary_linf}" \
        -v rel_tol="${GRADIENT_REL_TOL}" \
        -v abs_tol="${GRADIENT_ABS_TOL}" '
        function mag(value) {return value < 0 ? -value : value}
        function is_close(value, reference) {
            return mag(value - reference) <= abs_tol + rel_tol*mag(reference)
        }
        BEGIN {
            exit !(is_close(l2, ref_l2) && is_close(linf, ref_linf) &&
                is_close(boundary_l2, ref_boundary_l2) &&
                is_close(boundary_linf, ref_boundary_linf))
        }'
    then
        printf "PASS: serial %s: [%s] expected [%.8g %.8g %.8g %.8g]\n" \
            "${field}" "${row}" \
            "${ref_l2}" "${ref_linf}" "${ref_boundary_l2}" "${ref_boundary_linf}"
        return 0
    fi

    printf "FAIL: serial %s: [%s] expected [%.8g %.8g %.8g %.8g]\n" \
        "${field}" "${row}" \
        "${ref_l2}" "${ref_linf}" "${ref_boundary_l2}" "${ref_boundary_linf}"
    return 1
}

run_gradient_test() {
    if ! command -v Test-leastSquaresS4fGrad >/dev/null 2>&1; then
        echo "FAIL: Test-leastSquaresS4fGrad is not available"
        return 1
    fi

    prepare_case
    (
        cd "${CASE_DIR}"
        solids4Foam::convertCaseFormat . >/dev/null 2>&1
        solids4Foam::runApplication -o -s gradient blockMesh >/dev/null 2>&1
        sed -E -i.bak 's/^[[:space:]]*default[[:space:]]+none;/        default         leastSquaresS4f;/' system/fvSchemes
        rm -f system/fvSchemes.bak
        solids4Foam::runApplication -o -s serial Test-leastSquaresS4fGrad >/dev/null 2>&1
    )
}

run_high_order_grad_test() {
    local least_squares_type="$1"
    local suffix="${least_squares_type}.serial"
    local log_file="${CASE_DIR}/log.Test-highOrderGrad.${suffix}"

    if [[ -n "${FOAMEXTEND:-}" || "${WM_PROJECT_VERSION:-}" == "4.1" ]]; then
        echo "SKIP: Test-highOrderGrad is not available with foam-extend"
        return 0
    fi

    if ! command -v Test-highOrderGrad >/dev/null 2>&1; then
        echo "SKIP: Test-highOrderGrad is not available"
        return 0
    fi

    sed -E -i.bak \
        "s/^([[:space:]]*)type[[:space:]]+(movingLeastSquares|kExactLeastSquares);/\1type ${least_squares_type};/" \
        "${CASE_DIR}/constant/solidProperties.highOrder"
    rm -f "${CASE_DIR}/constant/solidProperties.highOrder.bak"

    if ! grep -qE \
        "^[[:space:]]*type[[:space:]]+${least_squares_type};" \
        "${CASE_DIR}/constant/solidProperties.highOrder"
    then
        echo "FAIL: Could not select ${least_squares_type}"
        return 1
    fi

    (
        cd "${CASE_DIR}"
        solids4Foam::runApplication \
            -o -s "${suffix}" Test-highOrderGrad >/dev/null 2>&1
    )

    if grep -q "Overall result: PASSED" "${log_file}" \
        && ! grep -q '^ERROR$' "${log_file}"
    then
        echo "PASS: Test-highOrderGrad (${least_squares_type}, serial)"
        return 0
    fi

    echo "FAIL: Test-highOrderGrad (${least_squares_type}, serial)"
    return 1
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
    rm -f "${CASE_DIR}/${GRADIENT_LOGFILE}"
    for approach in "${APPROACHES[@]}"; do
        echo
        echo "------------------------------------------------------------"
        echo "Testing approach: ${approach}"
        echo "------------------------------------------------------------"

        select_run_approach "${approach}"
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
        ( cd "${CASE_DIR}" && ./Allrun "${RUN_APPROACH}" > "${ALLRUN_LOGFILE}" 2>&1 )

        if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
            echo "Skipping ${approach} because it is unavailable in this environment"
            continue
        fi

        if [[ "${approach}" == highOrder-* ]]; then
            least_squares_type="${approach#highOrder-}"

            if ! run_high_order_grad_test "${least_squares_type}"; then
                failures=$((failures + 1))
            fi
        fi

        if ! check_top_right_disp "${approach}"; then
            failures=$((failures + 1))
        fi
    done

    echo
    echo "Testing leastSquaresS4f gradient calculation"
    if run_gradient_test; then
        if ! check_gradient_log; then
            failures=$((failures + 1))
        fi
    else
        failures=$((failures + 1))
    fi
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    if ! check_top_right_disp "check-only"; then
        failures=$((failures + 1))
    fi
    if ! check_gradient_log; then
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
