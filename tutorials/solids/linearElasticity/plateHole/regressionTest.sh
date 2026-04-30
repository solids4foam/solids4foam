#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# Plate-with-hole regression tests
# Checks numerical vs analytical solution for the displacement
# and pressure-displacement solution options.
# ============================================================

DISP_TOL=1e-7
POINT_DISP_TOL=1e-7
STRESS_TOL=2e5

PD_DISP_TOL=3.0e-4
PD_POINT_DISP_TOL=3.0e-4
PD_SIGMA_TOL=2.0e5
PD_P_TOL=9.0e4

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

PRESSURE_DISPLACEMENT_CASES=(
    "pressureDisplacementCompressible coarse"
    "pressureDisplacementCompressible medium"
    "pressureDisplacementIncompressible coarse"
    "pressureDisplacementIncompressible medium"
    "pressureDisplacementIncompressible poly"
)

echo "============================================================"
echo "Plate-with-hole regression tests"
echo "DDifference LInf        < ${DISP_TOL}"
echo "pointDDifference LInf   < ${POINT_DISP_TOL}"
echo "Stress component-0 LInf < ${STRESS_TOL}"
echo "Pressure-displacement DError max      < ${PD_DISP_TOL}"
echo "Pressure-displacement pointDError max < ${PD_POINT_DISP_TOL}"
echo "Pressure-displacement sigma*Err max   < ${PD_SIGMA_TOL}"
echo "Pressure-displacement pErr max        < ${PD_P_TOL}"
echo "============================================================"
echo

prepare_case() {
    local case_dir="$1"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        local base_item
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done
}

run_case() {
    local case_name="$1"
    shift

    local case_dir="${REGRESSION_ROOT}/${case_name}"

    prepare_case "${case_dir}"
    ( cd "${case_dir}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${case_dir}" && ./Allrun "$@" > "${ALLRUN_LOGFILE}" 2>&1 )

    echo "${case_dir}"
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
    rm -rf "${REGRESSION_ROOT}"
    mkdir -p "${REGRESSION_ROOT}"

    displacement_case_dir=$(run_case "displacement")

    for case_args in "${PRESSURE_DISPLACEMENT_CASES[@]}"; do
        IFS=' ' read -r approach mesh <<< "${case_args}"
        run_case "${approach}-${mesh}" "${approach}" "${mesh}" > /dev/null
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
    displacement_case_dir="${REGRESSION_ROOT}/displacement"
fi

extract_disp_linf() {
    local case_dir="$1"
    local field="$2"

    grep -A2 "Writing ${field} field" "${case_dir}/${SOLVER_LOGFILE}" \
        | grep "Norms:" -A1 \
        | tail -n 1 \
        | awk '{print $3}' \
        || true
}

extract_stress_linf_comp0() {
    local case_dir="$1"

    grep -A6 "Writing cellStressDifference field" \
        "${case_dir}/${SOLVER_LOGFILE}" \
        | tail -6 \
        | awk '
            /Component:[[:space:]]*0/ {getline; getline; print $3}
        ' \
        || true
}

extract_log_value() {
    local case_dir="$1"
    local label="$2"

    grep "${label}" "${case_dir}/${SOLVER_LOGFILE}" 2>/dev/null \
        | tail -n 1 \
        | awk '{print $NF}' \
        || true
}

check_less_than() {
    local case_name="$1"
    local label="$2"
    local value="$3"
    local tolerance="$4"

    if [[ -z "${value}" ]]; then
        echo "FAIL: ${case_name}: could not extract ${label}"
        failures=$((failures + 1))
    elif awk "BEGIN {exit !(${value} < ${tolerance})}"; then
        printf "PASS: %s: %s = %.6g\n" "${case_name}" "${label}" "${value}"
    else
        printf "FAIL: %s: %s = %.6g\n" "${case_name}" "${label}" "${value}"
        failures=$((failures + 1))
    fi
}

failures=0

if solids4Foam::regressionCaseSkipped \
    "${displacement_case_dir}/${ALLRUN_LOGFILE}"
then
    echo "SKIP: displacement"
else
    check_less_than \
        "displacement" "DDifference LInf" \
        "$(extract_disp_linf "${displacement_case_dir}" "DDifference")" \
        "${DISP_TOL}"
    check_less_than \
        "displacement" "pointDDifference LInf" \
        "$(extract_disp_linf "${displacement_case_dir}" "pointDDifference")" \
        "${POINT_DISP_TOL}"
    check_less_than \
        "displacement" "stress component-0 LInf" \
        "$(extract_stress_linf_comp0 "${displacement_case_dir}")" \
        "${STRESS_TOL}"
fi

for case_args in "${PRESSURE_DISPLACEMENT_CASES[@]}"; do
    IFS=' ' read -r approach mesh <<< "${case_args}"
    case_name="${approach}-${mesh}"
    case_dir="${REGRESSION_ROOT}/${case_name}"

    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: ${case_name}"
        continue
    fi

    check_less_than \
        "${case_name}" "DError, max" \
        "$(extract_log_value "${case_dir}" "DError, max")" "${PD_DISP_TOL}"
    check_less_than \
        "${case_name}" "pointDError, max" \
        "$(extract_log_value "${case_dir}" "pointDError, max")" \
        "${PD_POINT_DISP_TOL}"
    check_less_than \
        "${case_name}" "sigmaXXErr, max" \
        "$(extract_log_value "${case_dir}" "sigmaXXErr, max")" \
        "${PD_SIGMA_TOL}"
    check_less_than \
        "${case_name}" "sigmaXYErr, max" \
        "$(extract_log_value "${case_dir}" "sigmaXYErr, max")" \
        "${PD_SIGMA_TOL}"
    check_less_than \
        "${case_name}" "sigmaYYErr, max" \
        "$(extract_log_value "${case_dir}" "sigmaYYErr, max")" \
        "${PD_SIGMA_TOL}"
    check_less_than \
        "${case_name}" "pErr, max" \
        "$(extract_log_value "${case_dir}" "pErr, max")" "${PD_P_TOL}"
done

if [ "$CHECK_ONLY" = false ]; then
    for case_dir in "${REGRESSION_ROOT}"/*; do
        if [[ -d "${case_dir}" ]]; then
            ( cd "${case_dir}" && ./Allclean > /dev/null 2>&1 ) || true
        fi
    done
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
