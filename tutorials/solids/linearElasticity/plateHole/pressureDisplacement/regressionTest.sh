#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_ROOT="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# Plate-with-hole pressure-displacement regression tests
# Checks numerical vs analytical error summaries printed by setPlateHoleBC
# ============================================================

DISP_TOL=3.0e-4
POINT_DISP_TOL=3.0e-4
SIGMA_TOL=2.0e5
P_TOL=9.0e4

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

CASES=(
    "compressible/plateHoleHex/coarse"
    "compressible/plateHoleHex/medium"
    "incompressible/plateHoleHex/coarse"
    "incompressible/plateHoleHex/medium"
    "incompressible/plateHolePoly"
)

echo "============================================================"
echo "Plate-with-hole pressure-displacement regression tests"
echo "DError max       < ${DISP_TOL}"
echo "pointDError max  < ${POINT_DISP_TOL}"
echo "sigmaXXErr max   < ${SIGMA_TOL}"
echo "sigmaXYErr max   < ${SIGMA_TOL}"
echo "sigmaYYErr max   < ${SIGMA_TOL}"
echo "pErr max         < ${P_TOL}"
echo "============================================================"
echo

prepare_cases() {
    rm -rf "${CASE_ROOT}"
    mkdir -p "${CASE_ROOT}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_ROOT}/"
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
    prepare_cases

    for case_path in "${CASES[@]}"; do
        (
            cd "${CASE_ROOT}/${case_path}"
            ./Allclean > /dev/null 2>&1
        ) || true

        (
            cd "${CASE_ROOT}/${case_path}"
            ./Allrun > "${ALLRUN_LOGFILE}" 2>&1
        )
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

all_skipped=true

for case_path in "${CASES[@]}"; do
    if ! solids4Foam::regressionCaseSkipped \
        "${CASE_ROOT}/${case_path}/${ALLRUN_LOGFILE}"
    then
        all_skipped=false
    fi
done

if [ "${all_skipped}" = true ]; then
    echo "Skipping regression checks because all tutorials skipped"
    exit 0
fi

extract_log_value() {
    local case_path="$1"
    local label="$2"

    grep "${label}" "${CASE_ROOT}/${case_path}/${SOLVER_LOGFILE}" 2>/dev/null \
        | tail -n 1 \
        | awk '{print $NF}' \
        || true
}

check_less_than() {
    local case_path="$1"
    local label="$2"
    local value="$3"
    local tolerance="$4"

    if [[ -z "${value}" ]]; then
        echo "FAIL: ${case_path}: could not extract ${label}"
        failures=$((failures + 1))
    elif awk "BEGIN {exit !(${value} < ${tolerance})}"; then
        printf "PASS: %s: %s = %.6g\n" "${case_path}" "${label}" "${value}"
    else
        printf "FAIL: %s: %s = %.6g\n" "${case_path}" "${label}" "${value}"
        failures=$((failures + 1))
    fi
}

failures=0

for case_path in "${CASES[@]}"; do
    if solids4Foam::regressionCaseSkipped \
        "${CASE_ROOT}/${case_path}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: ${case_path}"
        continue
    fi

    check_less_than \
        "${case_path}" "DError, max" \
        "$(extract_log_value "${case_path}" "DError, max")" "${DISP_TOL}"
    check_less_than \
        "${case_path}" "pointDError, max" \
        "$(extract_log_value "${case_path}" "pointDError, max")" \
        "${POINT_DISP_TOL}"
    check_less_than \
        "${case_path}" "sigmaXXErr, max" \
        "$(extract_log_value "${case_path}" "sigmaXXErr, max")" \
        "${SIGMA_TOL}"
    check_less_than \
        "${case_path}" "sigmaXYErr, max" \
        "$(extract_log_value "${case_path}" "sigmaXYErr, max")" \
        "${SIGMA_TOL}"
    check_less_than \
        "${case_path}" "sigmaYYErr, max" \
        "$(extract_log_value "${case_path}" "sigmaYYErr, max")" \
        "${SIGMA_TOL}"
    check_less_than \
        "${case_path}" "pErr, max" \
        "$(extract_log_value "${case_path}" "pErr, max")" "${P_TOL}"
done

if [ "$CHECK_ONLY" = false ]; then
    for case_path in "${CASES[@]}"; do
        (
            cd "${CASE_ROOT}/${case_path}"
            ./Allclean > /dev/null 2>&1
        ) || true
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
