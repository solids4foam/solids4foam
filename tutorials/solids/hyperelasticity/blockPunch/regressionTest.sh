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
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"

APPROACHES=(
    segregated
    petscSnes
    highOrder
)

echo "============================================================"
echo "blockPunch regression test"
echo "Time steps: ${EXPECTED_TIME_STEPS}"
echo "Final time: ${REGRESSION_END_TIME}"
echo "Point-A disp_z between ${DISP_Z_MAX} and ${DISP_Z_MIN}"
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

    n_steps=$(awk '!/^#/ && NF {count++} END {print count + 0}' \
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

    if [ "$CHECK_ONLY" = false ]; then
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    fi
done

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
