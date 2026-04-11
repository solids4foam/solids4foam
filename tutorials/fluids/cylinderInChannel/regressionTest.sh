#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
PETSC_RETRY_CASE_DIR="${REGRESSION_ROOT}/petscRetry"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# cylinderInChannel regression test
# Uses the force history as a cheap fluid benchmark check.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

REG_END_TIME=5
DRAG_MIN=0.622
DRAG_MAX=0.624
LIFT_MIN=0.005
LIFT_MAX=0.012

ALLRUN_LOGFILE="log.Allrun"
PETSC_RETRY_LOGFILE="log.Allrun.petscRetry"

echo "============================================================"
echo "cylinderInChannel regression test"
echo "Regression end time = ${REG_END_TIME}"
echo "Final drag in [${DRAG_MIN}, ${DRAG_MAX}]"
echo "Final lift in [${LIFT_MIN}, ${LIFT_MAX}]"
echo "============================================================"
echo

prepare_case() {
    local dst_case="${1:-${CASE_DIR}}"

    rm -rf "${dst_case}"
    mkdir -p "${dst_case}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${dst_case}/"
    done

    sed -i.bak 's/^endTime[[:space:]]\+50;/endTime         5;/' "${dst_case}/system/controlDict"
    rm -f "${dst_case}/system/controlDict.bak"
}

find_force_file() {
    local candidate
    for candidate in \
        "${CASE_DIR}/forces/0/forces.dat" \
        "${CASE_DIR}/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/forces/0/forces.dat" \
        "${CASE_DIR}/postProcessing/fluid/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/fluid/forces/0/forces.dat"
    do
        if [[ -f "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done
    return 1
}

configure_petsc_retry_case() {
    local control_dict="${PETSC_RETRY_CASE_DIR}/system/controlDict"
    local fluid_properties="${PETSC_RETRY_CASE_DIR}/constant/fluidProperties.newtonIcoFluid"
    local fv_solution="${PETSC_RETRY_CASE_DIR}/system/fvSolution.newtonIcoFluid"

    awk '
    /^endTime[[:space:]]/ {
        print "endTime         0.11;"
        next
    }
    /^deltaT[[:space:]]/ {
        print "deltaT          0.1;"
        next
    }
    /^adjustTimeStep[[:space:]]/ {
        print "adjustTimeStep yes;"
        print ""
        print "maxDeltaT       0.1;"
        print ""
        print "minDeltaT       0.0001;"
        print ""
        print "logTimeStepAdjustments false;"
        next
    }
    { print }
    ' "${control_dict}" > "${control_dict}.tmp"
    mv "${control_dict}.tmp" "${control_dict}"

    awk '
    /^[[:space:]]*predictor[[:space:]]/ && !inserted {
        print
        print ""
        print "    stopOnPetscError false;"
        print "    maxTimeStepRetries 8;"
        inserted = 1
        next
    }
    { print }
    ' "${fluid_properties}" > "${fluid_properties}.tmp"
    mv "${fluid_properties}.tmp" "${fluid_properties}"

    awk '
    /^[[:space:]]*snes_max_it[[:space:]]/ {
        print "            snes_max_it \"3\";"
        next
    }
    { print }
    ' "${fv_solution}" > "${fv_solution}.tmp"
    mv "${fv_solution}.tmp" "${fv_solution}"
}

run_petsc_retry_check() {
    local retry_log="${PETSC_RETRY_CASE_DIR}/${PETSC_RETRY_LOGFILE}"
    local solver_log="${PETSC_RETRY_CASE_DIR}/log.solids4Foam"
    local run_status=0
    local retry_count=0
    local same_deltaT_retry_count=0
    local reduced_deltaT_retry_count=0
    local max_it_count=0
    local reset_count=0
    local failures=0

    echo "============================================================"
    echo "newtonIcoFluid PETSc retry regression check"
    echo "============================================================"

    if [[ -z "${PETSC_DIR:-}" ]]; then
        echo "Skipping this case as PETSc is not installed"
        echo "Please set the PETSC_DIR variable and rebuild solids4foam"
        return 0
    fi

    prepare_case "${PETSC_RETRY_CASE_DIR}"
    configure_petsc_retry_case

    ( cd "${PETSC_RETRY_CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

    set +e
    ( cd "${PETSC_RETRY_CASE_DIR}" && ./Allrun newtonIcoFluid > "${PETSC_RETRY_LOGFILE}" 2>&1 )
    run_status=$?
    set -e

    if solids4Foam::regressionCaseSkipped "${retry_log}"; then
        echo "Skipping PETSc retry check because the tutorial skipped in this environment"
        return 0
    fi

    if [[ ! -f "${solver_log}" ]]; then
        echo "FAIL: Could not find ${solver_log}"
        return 1
    fi

    if [[ "${run_status}" -ne 0 ]]; then
        printf "INFO: Allrun returned %d\n" "${run_status}"
    fi

    retry_count=$(grep -c "Retrying the failed PETSc time step" "${solver_log}" || true)
    same_deltaT_retry_count=$(grep -c "unchanged deltaT" "${solver_log}" || true)
    reduced_deltaT_retry_count=$(grep -c "with deltaT =" "${solver_log}" || true)
    max_it_count=$(grep -c "DIVERGED_MAX_IT" "${solver_log}" || true)
    reset_count=$(grep -c "Resetting PETSc SNES/KSP retry state" "${solver_log}" || true)

    if [[ "${retry_count}" -ge 2 ]]; then
        printf "PASS: observed %d PETSc retry attempts\n" "${retry_count}"
    else
        printf "FAIL: observed %d PETSc retry attempts\n" "${retry_count}"
        failures=$((failures + 1))
    fi

    if [[ "${same_deltaT_retry_count}" -ge 1 ]]; then
        printf "PASS: observed %d same-deltaT PETSc reset retries\n" \
            "${same_deltaT_retry_count}"
    else
        printf "FAIL: observed %d same-deltaT PETSc reset retries\n" \
            "${same_deltaT_retry_count}"
        failures=$((failures + 1))
    fi

    if [[ "${reduced_deltaT_retry_count}" -ge 1 ]]; then
        printf "PASS: observed %d reduced-deltaT PETSc retries\n" \
            "${reduced_deltaT_retry_count}"
    else
        printf "FAIL: observed %d reduced-deltaT PETSc retries\n" \
            "${reduced_deltaT_retry_count}"
        failures=$((failures + 1))
    fi

    if [[ "${max_it_count}" -ge 2 ]]; then
        printf "PASS: observed %d DIVERGED_MAX_IT reports\n" "${max_it_count}"
    else
        printf "FAIL: observed %d DIVERGED_MAX_IT reports\n" "${max_it_count}"
        failures=$((failures + 1))
    fi

    if [[ "${reset_count}" -ge 2 ]]; then
        printf "PASS: observed %d PETSc state resets\n" "${reset_count}"
    else
        printf "FAIL: observed %d PETSc state resets\n" "${reset_count}"
        failures=$((failures + 1))
    fi

    if grep -q "^End$" "${solver_log}"; then
        echo "PASS: retry path reached End"
    else
        echo "FAIL: retry path did not reach End"
        failures=$((failures + 1))
    fi

    if grep -Eq "Segmentation Violation|MPI_ABORT|Caught signal|^ERROR$" "${solver_log}"; then
        echo "FAIL: PETSc retry path crashed or returned an application error"
        failures=$((failures + 1))
    else
        echo "PASS: PETSc retry path did not crash"
    fi

    return "${failures}"
}

# ------------------------------------------------------------
# Clean & run case
# ------------------------------------------------------------

CHECK_ONLY=false
PETSC_RETRY_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run)
            CHECK_ONLY=true
            ;;
        --petsc-retry-only)
            PETSC_RETRY_ONLY=true
            ;;
        *)
            ;;
    esac
done

if [ "${PETSC_RETRY_ONLY}" = true ]; then
    run_petsc_retry_check
    exit $?
fi

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
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

force_file=""
if ! force_file=$(find_force_file); then
    echo "FAIL: Could not find force history output"
    exit 1
fi

final_drag=$(awk '
    END {
        gsub(/[()]/, "", $0)

        if (NF >= 13)
        {
            print $2 + $5
        }
        else if (NF >= 9)
        {
            print $2
        }
    }
' "${force_file}")
final_lift=$(awk '
    END {
        gsub(/[()]/, "", $0)

        if (NF >= 13)
        {
            print $3 + $6
        }
        else if (NF >= 9)
        {
            print $3
        }
    }
' "${force_file}")

if [[ -z "${final_drag}" || -z "${final_lift}" ]]; then
    echo "FAIL: Could not extract final drag/lift"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${final_drag} >= ${DRAG_MIN} && ${final_drag} <= ${DRAG_MAX})}"; then
    printf "PASS: Final drag = %.6g\n" "${final_drag}"
else
    printf "FAIL: Final drag = %.6g\n" "${final_drag}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_lift} >= ${LIFT_MIN} && ${final_lift} <= ${LIFT_MAX})}"; then
    printf "PASS: Final lift = %.6g\n" "${final_lift}"
else
    printf "FAIL: Final lift = %.6g\n" "${final_lift}"
    failures=$((failures + 1))
fi

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

    if ! run_petsc_retry_check; then
        failures=$((failures + 1))
    fi
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
