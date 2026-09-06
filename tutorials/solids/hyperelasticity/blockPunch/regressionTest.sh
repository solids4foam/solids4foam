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
# blockPunch regression test
# Checks four load steps with the segregated, PETSc SNES, and high-order
# approaches.
# ============================================================

REGRESSION_END_TIME=0.4
EXPECTED_TIME_STEPS=4
DISP_Z_MIN=-0.16
DISP_Z_MAX=-0.163

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
CONSTITUTIVE_LOGFILE="log.Test-mechanicalConstitutiveLaw"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"

APPROACHES=(
    segregated
    segregatedManager
    petscSnes
    highOrder
)

echo "============================================================"
echo "blockPunch regression test"
echo "Time steps: ${EXPECTED_TIME_STEPS}"
echo "Final time: ${REGRESSION_END_TIME}"
echo "Point-A disp_z between ${DISP_Z_MAX} and ${DISP_Z_MIN}"
echo "Plus the finite-strain mechanicalConstitutiveLaw framework checks"
echo "============================================================"
echo

prepare_case() {
    local approach="$1"
    local case_dir="${REGRESSION_ROOT}/${approach}"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done

    sed -i.bak \
        "s/^endTime[[:space:]][[:space:]]*[^;][^;]*;/endTime         ${REGRESSION_END_TIME};/" \
        "${case_dir}/system/controlDict"
    rm -f "${case_dir}/system/controlDict.bak"
}

extract_final_disp_z() {
    local case_dir="$1"

    awk 'END {print $4}' "${case_dir}/${DISP_FILE}" 2>/dev/null || true
}

check_history() {
    local case_dir="$1"
    local approach="$2"
    local n_steps
    local final_time

    # OpenFOAM.org writes an entry for the initial time whereas OpenFOAM.com
    # does not, so only the entries after time zero are counted
    n_steps=$(awk '!/^#/ && NF && $1 + 0 != 0 {count++} END {print count + 0}' \
        "${case_dir}/${DISP_FILE}" 2>/dev/null || true)
    final_time=$(awk '!/^#/ && NF {time=$1} END {print time}' \
        "${case_dir}/${DISP_FILE}" 2>/dev/null || true)

    if [[ "${n_steps}" != "${EXPECTED_TIME_STEPS}" ]]; then
        echo "FAIL: ${approach} wrote ${n_steps} displacement entries; expected ${EXPECTED_TIME_STEPS}"
        return 1
    fi

    if ! awk "BEGIN {exit !((${final_time} - ${REGRESSION_END_TIME})^2 < 1e-20)}"; then
        echo "FAIL: ${approach} final time is ${final_time}; expected ${REGRESSION_END_TIME}"
        return 1
    fi

    return 0
}

check_log() {
    local case_dir="$1"
    local approach="$2"

    if grep -Eq "DIVERGED|did not converge|FOAM FATAL ERROR" \
        "${case_dir}/${SOLVER_LOGFILE}"; then
        echo "FAIL: Solver log reports non-convergence for ${approach}"
        return 1
    fi

    return 0
}

check_displacement() {
    local case_dir="$1"
    local approach="$2"
    local disp_z

    disp_z=$(extract_final_disp_z "${case_dir}")

    if [[ -z "${disp_z}" ]]; then
        echo "FAIL: Could not extract point-A displacement for ${approach}"
        return 1
    fi

    if awk "BEGIN {exit !(${disp_z} <= ${DISP_Z_MIN} && ${disp_z} >= ${DISP_Z_MAX})}"; then
        printf "PASS: %s final point-A disp_z = %.7g\n" \
            "${approach}" "${disp_z}"
        return 0
    else
        printf "FAIL: %s final point-A disp_z = %.7g; expected between %.7g and %.7g\n" \
            "${approach}" "${disp_z}" "${DISP_Z_MAX}" "${DISP_Z_MIN}"
        return 1
    fi
}

# Exercise the finite-strain paths of the mechanicalConstitutiveLawManager on
# this case. This tutorial is used because its law is neoHookeanElastic, which
# implements no small-strain evaluation, so it is the only runtime coverage of
# the finite-strain kinematics, the finite-difference fourth-order tangent, and
# the boundary integration points of a face-centred topology
run_constitutive_test() {
    local case_dir="$1"

    if ! command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
        echo "SKIP: Test-mechanicalConstitutiveLaw not found in PATH"
        return 0
    fi

    if [[ ! -d "${case_dir}/constant/polyMesh" ]]; then
        echo "SKIP: mechanicalConstitutiveLaw checks (case has no mesh)"
        return 0
    fi

    if ( cd "${case_dir}" && Test-mechanicalConstitutiveLaw \
            > "${CONSTITUTIVE_LOGFILE}" 2>&1 )
    then
        local n_passed
        n_passed=$(grep -c 'PASS:' "${case_dir}/${CONSTITUTIVE_LOGFILE}" || true)

        if (( n_passed == 0 )); then
            echo "SKIP: mechanicalConstitutiveLaw checks (no checks reported)"
            return 0
        fi

        echo "PASS: mechanicalConstitutiveLaw checks (${n_passed} checks)"
        return 0
    fi

    echo "FAIL: mechanicalConstitutiveLaw checks"
    grep 'FAIL:' "${case_dir}/${CONSTITUTIVE_LOGFILE}" || true
    return 1
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
    if ! command -v gmsh > /dev/null 2>&1; then
        echo "Skipping regression checks because Gmsh is not installed"
        exit 0
    fi
fi

failures=0
declare -A RESULT_DISP

# The framework checks depend only on the mesh and the material, so they are
# run once rather than once per approach
constitutive_tested=false

for approach in "${APPROACHES[@]}"; do
    CASE_DIR="${REGRESSION_ROOT}/${approach}"

    echo
    echo "------------------------------------------------------------"
    echo "Testing approach: ${approach}"
    echo "------------------------------------------------------------"

    if [ "$CHECK_ONLY" = false ]; then
        prepare_case "${approach}"
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
        ( cd "${CASE_DIR}" && ./Allrun "${approach}" > "${ALLRUN_LOGFILE}" 2>&1 )
    else
        echo "Running in check-only mode: skipping Allclean and Allrun"
    fi

    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping ${approach} because it is unavailable in this environment"
        continue
    fi

    if ! check_log "${CASE_DIR}" "${approach}"; then
        failures=$((failures + 1))
        continue
    fi

    if ! check_history "${CASE_DIR}" "${approach}"; then
        failures=$((failures + 1))
        continue
    fi

    if ! check_displacement "${CASE_DIR}" "${approach}"; then
        failures=$((failures + 1))
    fi

    RESULT_DISP["${approach}"]=$(extract_final_disp_z "${CASE_DIR}")

    # Before the Allclean below, which removes the mesh
    if [ "${constitutive_tested}" = false ]; then
        constitutive_tested=true

        if ! run_constitutive_test "${CASE_DIR}"; then
            failures=$((failures + 1))
        fi
    fi

    if [ "$CHECK_ONLY" = false ]; then
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    fi
done

# segregatedManager is solidProperties.segregated plus the framework switch and
# nothing else, so the two must agree. This is the only end-to-end comparison
# of the framework's neoHookeanElastic against the legacy law
if [[ -n "${RESULT_DISP[segregated]:-}" \
   && -n "${RESULT_DISP[segregatedManager]:-}" ]]
then
    a="${RESULT_DISP[segregated]}"
    b="${RESULT_DISP[segregatedManager]}"

    if awk "BEGIN {exit !(($a - $b)^2 <= (1e-6*$a)^2)}"; then
        printf "PASS: legacy and framework disp_z agree (%.8g vs %.8g)\n" \
            "$a" "$b"
    else
        printf "FAIL: legacy and framework disp_z differ (%.8g vs %.8g)\n" \
            "$a" "$b"
        failures=$((failures + 1))
    fi
else
    echo "SKIP: cross-check needs both segregated approaches to have run"
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
