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
# Hyperelastic cooksMembrane regression test
# Checks the vertical displacement of the reference point, the mid-point of
# the loaded edge, at the end of the loading.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

# solids4foam gives 0.0146148 m on the 24 x 24 mesh with both approaches
# (OpenFOAM-v2512). The band is +/- 0.5% about this value. For comparison,
# Abaqus CPE4H gives 0.0146183 m on a mesh with the same number of cells
# (576), and the mesh-converged deal.II Q2 value is 0.01474 m; see
# reference/abaqus_tipDisplacement.dat and reference/dealII_tipDisplacement.dat
REF_DISP_MIN=0.01454
REF_DISP_MAX=0.01469

# The load is fully applied at the end time
END_TIME=30

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"

APPROACHES=(
    segregated
    petscSnes
)

echo "============================================================"
echo "cooksMembrane (hyperelastic) regression test"
echo "Reference point displacement in [${REF_DISP_MIN}, ${REF_DISP_MAX}] m"
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

check_solver_log() {
    local approach="$1"
    local log_file="${CASE_DIR}/${SOLVER_LOGFILE}"

    if [[ ! -f "${log_file}" ]]; then
        echo "FAIL: Could not find ${SOLVER_LOGFILE} for ${approach}"
        return 1
    fi

    if grep -qE 'FOAM FATAL|^ERROR$|\[stack trace\]' "${log_file}" \
        || ! grep -q '^End' "${log_file}"
    then
        echo "FAIL: solids4Foam did not complete for ${approach}"
        return 1
    fi

    return 0
}

check_ref_disp() {
    local approach="$1"
    local value_file
    local final_time
    local ref_disp

    value_file=$(find "${CASE_DIR}/postProcessing" -name 'solidPointDisplacement_pointDisp.dat' -print 2>/dev/null | tail -n 1)
    if [[ -z "${value_file}" ]]; then
        echo "FAIL: Could not find point displacement output for ${approach}"
        return 1
    fi

    final_time=$(awk 'END {print $1}' "${value_file}")
    ref_disp=$(awk 'END {print $3}' "${value_file}")

    if [[ -z "${ref_disp}" ]]; then
        echo "FAIL: Could not extract reference point displacement for ${approach}"
        return 1
    fi

    if ! awk "BEGIN {exit !(${final_time} == ${END_TIME})}"; then
        echo "FAIL: Final time = ${final_time}, expected ${END_TIME} (${approach})"
        return 1
    fi

    if awk "BEGIN {exit !(${ref_disp} >= ${REF_DISP_MIN} && ${ref_disp} <= ${REF_DISP_MAX})}"; then
        printf "PASS: Reference point displacement = %.6g (%s)\n" "${ref_disp}" "${approach}"
        return 0
    fi

    printf "FAIL: Reference point displacement = %.6g (%s)\n" "${ref_disp}" "${approach}"
    return 1
}

failures=0

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
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

        if ! check_solver_log "${approach}"; then
            failures=$((failures + 1))
            continue
        fi

        if ! check_ref_disp "${approach}"; then
            failures=$((failures + 1))
        fi
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    if ! check_solver_log "check-only"; then
        failures=$((failures + 1))
    elif ! check_ref_disp "check-only"; then
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
