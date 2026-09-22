#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
FRAMEWORK_DIR="${REGRESSION_ROOT}/framework"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# thermalCavity regression test
# Uses the fsiConvergenceData output as a cheap convergence check.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

N_FSI_CORRECTORS_MAX=50

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "thermalCavity regression test"
echo "nFsiCorrectors < ${N_FSI_CORRECTORS_MAX}"
echo "============================================================"
echo

prepare_case() {
    local d="${1:-${CASE_DIR}}"
    rm -rf "${d}"
    mkdir -p "${d}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${d}/"
    done
}

shorten_case() {
    local controlDict="${1:-${CASE_DIR}}/system/controlDict"
    sed -i.bak 's/^endTime[[:space:]]\+10;/endTime         0.1;/' "${controlDict}"
    rm -f "${controlDict}.bak"
}

find_fsi_data() {
    local root="${1:-${CASE_DIR}}"
    local candidate
    for candidate in \
        "${root}/postProcessing/0/fsiConvergenceData.dat" \
        "${root}/postProcessing/fluid/0/fsiConvergenceData.dat" \
        "${root}/postProcessing/solid/0/fsiConvergenceData.dat"
    do
        if [[ -f "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done

    find "${root}/postProcessing" -name 'fsiConvergenceData.dat' \
        -print 2>/dev/null | tail -n 1
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
    prepare_case "${CASE_DIR}"
    shorten_case "${CASE_DIR}"
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )

    # The framework arm differs in this one entry and nothing else. thermalSolid
    # asks the constitutive implementation only for the density, so this checks
    # that the framework supplies the same rho the legacy model did
    prepare_case "${FRAMEWORK_DIR}"
    shorten_case "${FRAMEWORK_DIR}"
    sed -i.bak \
        's|^\( *\)solutionTolerance|\1useMechanicalConstitutiveLawManager yes;\n\1solutionTolerance|' \
        "${FRAMEWORK_DIR}/constant/solid/solidProperties"
    rm -f "${FRAMEWORK_DIR}/constant/solid/solidProperties.bak"

    if ! grep -q "useMechanicalConstitutiveLawManager" \
        "${FRAMEWORK_DIR}/constant/solid/solidProperties"
    then
        echo "FAIL: could not set the framework switch on the framework arm"
        exit 1
    fi

    ( cd "${FRAMEWORK_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

fsi_data=$(find_fsi_data)
if [[ -z "${fsi_data}" ]]; then
    echo "FAIL: Could not find fsiConvergenceData output"
    exit 1
fi

n_fsi_correctors=$(grep -v '^[[:space:]]*#' "${fsi_data}" | tail -n 1 | awk '{print $2}')

if [[ -z "${n_fsi_correctors}" ]]; then
    echo "FAIL: Could not extract nFsiCorrectors"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${n_fsi_correctors} < ${N_FSI_CORRECTORS_MAX})}"; then
    printf "PASS: nFsiCorrectors = %.6g\n" "${n_fsi_correctors}"
else
    printf "FAIL: nFsiCorrectors = %.6g\n" "${n_fsi_correctors}"
    failures=$((failures + 1))
fi

# ------------------------------------------------------------
# The framework arm
#
# thermalSolid takes only the density from the constitutive implementation, so
# the two arms should agree exactly: the same rho drives the same conjugate
# heat transfer and the same coupling iteration count
# ------------------------------------------------------------

if [ "$CHECK_ONLY" = false ]; then
    if solids4Foam::regressionCaseSkipped "${FRAMEWORK_DIR}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: framework arm skipped in this environment"
    else
        fw_log=$(find "${FRAMEWORK_DIR}" -name "${SOLVER_LOGFILE}" | tail -n 1)
        main_log=$(find "${CASE_DIR}" -name "${SOLVER_LOGFILE}" | tail -n 1)

        # Each arm must have taken the path it was set up for, or this
        # compares the legacy path against itself and proves nothing
        if [[ -n "${main_log}" ]] \
            && grep -q "mechanicalConstitutiveLawManager" "${main_log}"
        then
            echo "FAIL: the legacy arm used the framework"
            failures=$((failures + 1))
        else
            echo "PASS: legacy arm took the legacy path"
        fi

        if [[ -n "${fw_log}" ]] \
            && grep -q "mechanicalConstitutiveLawManager" "${fw_log}"
        then
            echo "PASS: framework arm took the framework path"
        else
            echo "FAIL: framework arm did not take the framework path"
            failures=$((failures + 1))
        fi

        fw_data=$(find_fsi_data "${FRAMEWORK_DIR}")
        if [[ -z "${fw_data}" ]]; then
            echo "FAIL: framework arm produced no fsiConvergenceData"
            failures=$((failures + 1))
        elif [[ "${fw_data}" == "${fsi_data}" ]]; then
            # find_fsi_data once had CASE_DIR hard-coded in its preferred
            # candidates and ignored its argument, so this returned the legacy
            # file and the comparison below was legacy against itself. It
            # passed. Assert the two arms are read from two different files
            echo "FAIL: both arms read the same fsiConvergenceData file,"
            echo "      so the comparison below would be vacuous:"
            echo "      ${fw_data}"
            failures=$((failures + 1))
        else
            fw_n=$(grep -v '^[[:space:]]*#' "${fw_data}" | tail -n 1 \
                | awk '{print $2}')

            if [[ "${fw_n}" == "${n_fsi_correctors}" ]]; then
                printf "PASS: framework nFsiCorrectors = %s, as legacy\n" \
                    "${fw_n}"
            else
                printf "FAIL: nFsiCorrectors differ (%s legacy, %s framework)\n" \
                    "${n_fsi_correctors}" "${fw_n}"
                failures=$((failures + 1))
            fi
        fi
    fi
fi

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${FRAMEWORK_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
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
