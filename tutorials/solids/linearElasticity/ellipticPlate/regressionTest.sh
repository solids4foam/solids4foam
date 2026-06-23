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

EPS_MIN=5.84e-5
EPS_MAX=5.86e-5
SIGMA_MIN=1.41e7
SIGMA_MAX=1.43e7

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

APPROACHES=(
    segregated
    petscSnes
)

# Settings for Test-leastSquaresS4fGrad test
GRADIENT_REL_TOL=1e-5
GRADIENT_ABS_TOL=1e-10
GRADIENT_SERIAL_LOGFILE="log.Test-leastSquaresS4fGrad.serial"
GRADIENT_PARALLEL_LOGFILE="log.Test-leastSquaresS4fGrad.parallel"

echo "============================================================"
echo "ellipticPlate regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
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

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    rm -f "${CASE_DIR}/${GRADIENT_SERIAL_LOGFILE}" "${CASE_DIR}/${GRADIENT_PARALLEL_LOGFILE}"
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

check_gradient_log() {
    local logfile="$1"
    local mode="$2"
    local failures=0

    check_gradient_reference_field "${logfile}" "${mode}" phi_x 0 0 0 0 \
        || failures=$((failures + 1))
    check_gradient_reference_field "${logfile}" "${mode}" phi_linear 0 0 0 0 \
        || failures=$((failures + 1))
    check_gradient_reference_field "${logfile}" "${mode}" \
        phi_quadratic 0.0449495 0.302343 0.078008 0.302343 \
        || failures=$((failures + 1))
    check_gradient_reference_field "${logfile}" "${mode}" \
        phi_quadratic_cross 0.00882968 0.0421144 0.0142969 0.0421144 \
        || failures=$((failures + 1))
    check_gradient_reference_field "${logfile}" "${mode}" \
        phi_trigonometric 1.56111 5.49017 1.86558 5.49017 \
        || failures=$((failures + 1))

    return "${failures}"
}

extract_gradient_row() {
    local logfile="$1"
    local field="$2"

    awk -v field="${field}" \
        '$1 == field && NF >= 5 {print $2, $3, $4, $5; exit}' \
        "${logfile}" 2>/dev/null || true
}

check_gradient_reference_field() {
    local logfile="$1"
    local mode="$2"
    local field="$3"
    local ref_l2="$4"
    local ref_linf="$5"
    local ref_boundary_l2="$6"
    local ref_boundary_linf="$7"
    local row
    local l2 linf boundary_l2 boundary_linf

    row=$(extract_gradient_row "${logfile}" "${field}")
    if [[ -z "${row}" ]]; then
        echo "FAIL: Could not extract ${mode} gradient norms for ${field}"
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
        printf "PASS: %s %s: [%s] expected [%.8g %.8g %.8g %.8g]\n" \
            "${mode}" "${field}" "${row}" \
            "${ref_l2}" "${ref_linf}" "${ref_boundary_l2}" "${ref_boundary_linf}"
        return 0
    fi

    printf "FAIL: %s %s: [%s] expected [%.8g %.8g %.8g %.8g]\n" \
        "${mode}" "${field}" "${row}" \
        "${ref_l2}" "${ref_linf}" "${ref_boundary_l2}" "${ref_boundary_linf}"
    return 1
}

run_gradient_tests() {
    if ! command -v Test-leastSquaresS4fGrad >/dev/null 2>&1; then
        echo "FAIL: Test-leastSquaresS4fGrad is not available"
        return 1
    fi

    prepare_case

    (
        cd "${CASE_DIR}"
        solids4Foam::runApplication -o -s gradient blockMesh >/dev/null 2>&1
        sed -E -i.bak 's/^[[:space:]]*default[[:space:]]+none;/        default         leastSquaresS4f;/' system/fvSchemes
        rm -f system/fvSchemes.bak
        solids4Foam::runApplication -o -s serial Test-leastSquaresS4fGrad >/dev/null 2>&1

        solids4Foam::runApplication -o -s gradient decomposePar -force >/dev/null 2>&1
        solids4Foam::runParallel -o -s parallel -n 2 \
            Test-leastSquaresS4fGrad >/dev/null 2>&1
    )
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

    echo
    echo "Testing leastSquaresS4f gradients in serial and parallel"
    if run_gradient_tests; then
        if ! check_gradient_log "${CASE_DIR}/${GRADIENT_SERIAL_LOGFILE}" serial; then
            failures=$((failures + 1))
        fi
        if ! check_gradient_log "${CASE_DIR}/${GRADIENT_PARALLEL_LOGFILE}" parallel; then
            failures=$((failures + 1))
        fi
    else
        failures=$((failures + 1))
    fi
else
    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    if ! check_solver_extrema "check-only"; then
        failures=$((failures + 1))
    fi
    if ! check_gradient_log "${CASE_DIR}/${GRADIENT_SERIAL_LOGFILE}" serial; then
        failures=$((failures + 1))
    fi
    if ! check_gradient_log "${CASE_DIR}/${GRADIENT_PARALLEL_LOGFILE}" parallel; then
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
