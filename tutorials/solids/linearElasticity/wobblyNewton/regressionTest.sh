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
# wobblyNewton regression test
# Checks Test-leastSquaresS4fGrad on the default mesh.
# ============================================================

# Settings for Test-leastSquaresS4fGrad test
GRADIENT_MAX_ERROR=1e-12
# Keep this at two ranks for CI runners with two available MPI slots.
GRADIENT_N_PROCS=2
GRADIENT_SERIAL_LOGFILE="log.Test-leastSquaresS4fGrad.serial"
GRADIENT_PARALLEL_LOGFILE="log.Test-leastSquaresS4fGrad.parallel"

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

if [ "$CHECK_ONLY" = false ] && ! command -v cartesianMesh >/dev/null 2>&1; then
    exit 0
fi

echo "============================================================"
echo "wobblyNewton regression test"
echo "Gradient maximum error: ${GRADIENT_MAX_ERROR}"
echo "Gradient parallel subdomains: ${GRADIENT_N_PROCS}"
echo "============================================================"
echo

extract_gradient_row() {
    local logfile="$1"
    local field="$2"

    awk -v field="${field}" \
        '$1 == field && NF >= 5 {print $2, $3, $4, $5; exit}' \
        "${logfile}" 2>/dev/null || true
}

check_gradient_log() {
    local logfile="$1"
    local mode="$2"

    if [[ ! -s "${logfile}" ]]; then
        echo "FAIL: ${mode} gradient log not found: ${logfile}"
        return 1
    fi

    if ! grep -q '^phi_' "${logfile}"; then
        echo "FAIL: ${mode} gradient result table not found in ${logfile}"
        return 1
    fi

    return 0
}

check_gradient_field_error() {
    local logfile="$1"
    local mode="$2"
    local field="$3"
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
        -v max_error="${GRADIENT_MAX_ERROR}" '
        function mag(value) {return value < 0 ? -value : value}
        BEGIN {
            exit !(mag(l2) <= max_error && mag(linf) <= max_error &&
                mag(boundary_l2) <= max_error &&
                mag(boundary_linf) <= max_error)
        }'
    then
        printf "PASS: %s %s: [%s]\n" "${mode}" "${field}" "${row}"
        return 0
    fi

    printf "FAIL: %s %s: [%s]\n" "${mode}" "${field}" "${row}"
    return 1
}

check_gradient_logs() {
    local failures=0

    check_gradient_log "${CASE_DIR}/${GRADIENT_SERIAL_LOGFILE}" serial \
        || failures=$((failures + 1))
    check_gradient_log "${CASE_DIR}/${GRADIENT_PARALLEL_LOGFILE}" parallel \
        || failures=$((failures + 1))

    if (( failures != 0 )); then
        return "${failures}"
    fi

    check_gradient_field_error "${CASE_DIR}/${GRADIENT_SERIAL_LOGFILE}" \
        serial phi_x || failures=$((failures + 1))
    check_gradient_field_error "${CASE_DIR}/${GRADIENT_SERIAL_LOGFILE}" \
        serial phi_linear || failures=$((failures + 1))
    check_gradient_field_error "${CASE_DIR}/${GRADIENT_PARALLEL_LOGFILE}" \
        parallel phi_x || failures=$((failures + 1))
    check_gradient_field_error "${CASE_DIR}/${GRADIENT_PARALLEL_LOGFILE}" \
        parallel phi_linear || failures=$((failures + 1))

    return "${failures}"
}

run_gradient_tests() {
    prepare_case
    rm -f "${CASE_DIR}/${GRADIENT_SERIAL_LOGFILE}" \
        "${CASE_DIR}/${GRADIENT_PARALLEL_LOGFILE}"

    if ! command -v cartesianMesh >/dev/null 2>&1; then
        echo "Skipping gradient test because cartesianMesh is not available"
        return 0
    fi

    if ! command -v Test-leastSquaresS4fGrad >/dev/null 2>&1; then
        echo "FAIL: Test-leastSquaresS4fGrad is not available"
        return 1
    fi

    (
        failures=0

        cd "${CASE_DIR}"
        solids4Foam::convertCaseFormat . >/dev/null 2>&1

        if ! solids4Foam::runApplication -o -s gradient \
            cartesianMesh >/dev/null 2>&1
        then
            echo "FAIL: cartesianMesh failed for leastSquaresS4f gradient test"
            failures=$((failures + 1))
        fi

        if ! solids4Foam::runApplication -o -s gradient \
            surfaceToPatch base.stl -noFunctionObjects >/dev/null 2>&1
        then
            echo "FAIL: surfaceToPatch failed for leastSquaresS4f gradient test"
            failures=$((failures + 1))
        fi

        if [[ -d 0.01/polyMesh ]]; then
            rm -rf constant/polyMesh
            mv 0.01/polyMesh constant/
            rm -rf 0.01
        else
            echo "FAIL: surfaceToPatch did not create 0.01/polyMesh"
            failures=$((failures + 1))
        fi

        if ! sed -E -i.bak \
            's/^[[:space:]]*default[[:space:]]+none;/        default         leastSquaresS4f;/' \
            system/fvSchemes
        then
            echo "FAIL: Could not enable leastSquaresS4f in system/fvSchemes"
            failures=$((failures + 1))
        fi
        rm -f system/fvSchemes.bak

        if ! sed -E -i.bak \
            "s/^[[:space:]]*numberOfSubdomains[[:space:]]+[0-9]+;/numberOfSubdomains ${GRADIENT_N_PROCS};/" \
            system/decomposeParDict
        then
            echo "FAIL: Could not set gradient test subdomain count"
            failures=$((failures + 1))
        fi
        rm -f system/decomposeParDict.bak

        if ! solids4Foam::runApplication -o -s serial \
            Test-leastSquaresS4fGrad >/dev/null 2>&1
        then
            echo "FAIL: serial Test-leastSquaresS4fGrad run failed"
            failures=$((failures + 1))
        fi

        if ! solids4Foam::runApplication -o -s gradient \
            decomposePar -force >/dev/null 2>&1
        then
            echo "FAIL: decomposePar failed for parallel leastSquaresS4f gradient test"
            failures=$((failures + 1))
        fi

        if ! (
            set +u
            solids4Foam::runParallel -o -s parallel -n "${GRADIENT_N_PROCS}" \
                Test-leastSquaresS4fGrad
        ) >/dev/null 2>&1
        then
            echo "FAIL: parallel Test-leastSquaresS4fGrad run failed"
            failures=$((failures + 1))
        fi

        exit "${failures}"
    )
}

failures=0

if [ "$CHECK_ONLY" = false ]; then
    echo "Testing leastSquaresS4f gradients in serial and parallel"
    if ! run_gradient_tests; then
        echo "One or more leastSquaresS4f gradient commands failed; checking produced logs"
    fi
    if ! check_gradient_logs; then
        failures=$((failures + 1))
    fi
else
    echo "Running in check-only mode: skipping mesh generation and gradient runs"
    if ! check_gradient_logs; then
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
