#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
PETSC_RETRY_CASE_DIR="${REGRESSION_ROOT}/petscRetry"
PETSC_RETRY_RECOVERY_CASE_DIR="${REGRESSION_ROOT}/petscRetryRecovery"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# cavityFlexibleBottom FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

DISP_MAX_TOL=0.005      # max displacement absolute tolerance
FORCE_MEAN_TOL=0.150    # mean force tolerance

# Number of samples from end of force.dat to average
FORCE_AVG_SAMPLES=40

# Reference values
REF_MAX_DISP=-0.21
REF_MEAN_FORCE=-5.150

# Log files
ALLRUN_LOGFILE="log.Allrun"
PETSC_RETRY_LOGFILE="log.Allrun.petscRetry"
PETSC_RETRY_RECOVERY_LOGFILE="log.Allrun.petscRetryRecovery"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

echo "============================================================"
echo "Beam-in-cross-flow FSI regression test"
echo "Max displacement difference < ${DISP_MAX_TOL}"
echo "Mean force difference       < ${FORCE_MEAN_TOL}"
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
}

configure_petsc_retry_case() {
    local control_dict="${PETSC_RETRY_CASE_DIR}/system/controlDict"
    local fsi_properties="${PETSC_RETRY_CASE_DIR}/constant/fsiProperties.monolithic"
    local fv_solution="${PETSC_RETRY_CASE_DIR}/system/fvSolution.monolithic"

    awk '
    /^endTime[[:space:]]/ {
        print "endTime         0.5;"
        next
    }
    /^adjustTimeStep[[:space:]]/ {
        print "adjustTimeStep  yes;"
        print ""
        print "maxDeltaT       0.5;"
        print ""
        print "minDeltaT       0.01;"
        print ""
        print "logTimeStepAdjustments false;"
        next
    }
    { print }
    ' "${control_dict}" > "${control_dict}.tmp"
    mv "${control_dict}.tmp" "${control_dict}"

    awk '
    /^[[:space:]]*passViscousStress[[:space:]]/ && !inserted {
        print "    stopOnPetscError false;"
        print "    maxTimeStepRetries 2;"
        print ""
        inserted = 1
    }
    { print }
    ' "${fsi_properties}" > "${fsi_properties}.tmp"
    mv "${fsi_properties}.tmp" "${fsi_properties}"

    awk '
    /^[[:space:]]*snes_stol[[:space:]]/ && !inserted {
        print
        print "            snes_max_it \"1\";"
        inserted = 1
        next
    }
    { print }
    ' "${fv_solution}" > "${fv_solution}.tmp"
    mv "${fv_solution}.tmp" "${fv_solution}"
}

configure_petsc_retry_recovery_case() {
    local control_dict="${PETSC_RETRY_RECOVERY_CASE_DIR}/system/controlDict"
    local fsi_properties="${PETSC_RETRY_RECOVERY_CASE_DIR}/constant/fsiProperties.monolithic"
    local fv_solution="${PETSC_RETRY_RECOVERY_CASE_DIR}/system/fvSolution.monolithic"

    awk '
    /^endTime[[:space:]]/ {
        print "endTime         5.1;"
        next
    }
    /^deltaT[[:space:]]/ {
        print "deltaT          5;"
        next
    }
    /^adjustTimeStep[[:space:]]/ {
        print "adjustTimeStep  yes;"
        print ""
        print "maxDeltaT       5;"
        print ""
        print "minDeltaT       0.001;"
        print ""
        print "logTimeStepAdjustments false;"
        next
    }
    { print }
    ' "${control_dict}" > "${control_dict}.tmp"
    mv "${control_dict}.tmp" "${control_dict}"

    awk '
    /^[[:space:]]*passViscousStress[[:space:]]/ && !inserted {
        print "    stopOnPetscError false;"
        print "    maxTimeStepRetries 8;"
        print ""
        inserted = 1
    }
    { print }
    ' "${fsi_properties}" > "${fsi_properties}.tmp"
    mv "${fsi_properties}.tmp" "${fsi_properties}"

    awk '
    /^[[:space:]]*snes_stol[[:space:]]/ && !inserted {
        print
        print "            snes_max_it \"2\";"
        inserted = 1
        next
    }
    { print }
    ' "${fv_solution}" > "${fv_solution}.tmp"
    mv "${fv_solution}.tmp" "${fv_solution}"
}

run_petsc_retry_check() {
    local run_status=0
    local retry_log="${PETSC_RETRY_CASE_DIR}/${PETSC_RETRY_LOGFILE}"
    local solver_log="${PETSC_RETRY_CASE_DIR}/log.solids4Foam"
    local retry_count=0
    local same_deltaT_retry_count=0
    local reduced_deltaT_retry_count=0
    local max_it_count=0
    local reset_count=0
    local failures=0

    echo "============================================================"
    echo "PETSc retry regression check"
    echo "============================================================"

    if [[ -z "${PETSC_DIR:-}" ]]; then
        echo "Skipping PETSc retry check because PETSc is not installed"
        echo "Please set the PETSC_DIR variable and rebuild solids4foam"
        return 0
    fi

    prepare_case "${PETSC_RETRY_CASE_DIR}"
    configure_petsc_retry_case

    ( cd "${PETSC_RETRY_CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

    set +e
    ( cd "${PETSC_RETRY_CASE_DIR}" && ./Allrun monolithic > "${PETSC_RETRY_LOGFILE}" 2>&1 )
    run_status=$?
    set -e

    if declare -F "solids4Foam::regressionCaseSkipped" > /dev/null \
        && solids4Foam::regressionCaseSkipped "${retry_log}"
    then
        echo "Skipping PETSc retry check because the tutorial skipped in this environment"
        return 0
    fi

    if [[ ! -f "${solver_log}" ]]; then
        echo "FAIL: Could not find ${solver_log}"
        return 1
    fi

    retry_count=$(grep -c "Retrying the failed PETSc time step" "${solver_log}" || true)
    same_deltaT_retry_count=$(grep -c "unchanged deltaT" "${solver_log}" || true)
    reduced_deltaT_retry_count=$(grep -c "with deltaT =" "${solver_log}" || true)
    max_it_count=$(grep -c "DIVERGED_MAX_IT" "${solver_log}" || true)
    reset_count=$(grep -c "Resetting PETSc SNES/KSP retry state" "${solver_log}" || true)

    if [[ "${retry_count}" -ge 3 ]]; then
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

    if [[ "${reduced_deltaT_retry_count}" -ge 2 ]]; then
        printf "PASS: observed %d reduced-deltaT PETSc retries\n" \
            "${reduced_deltaT_retry_count}"
    else
        printf "FAIL: observed %d reduced-deltaT PETSc retries\n" \
            "${reduced_deltaT_retry_count}"
        failures=$((failures + 1))
    fi

    if [[ "${max_it_count}" -ge 4 ]]; then
        printf "PASS: observed %d DIVERGED_MAX_IT reports\n" "${max_it_count}"
    else
        printf "FAIL: observed %d DIVERGED_MAX_IT reports\n" "${max_it_count}"
        failures=$((failures + 1))
    fi

    if grep -q "Exceeded the maximum number of failed PETSc retries (2)" "${solver_log}"; then
        echo "PASS: retry cap failure was reported"
    else
        echo "FAIL: retry cap failure was not reported"
        failures=$((failures + 1))
    fi

    if [[ "${reset_count}" -ge 3 ]]; then
        printf "PASS: observed %d PETSc state resets\n" "${reset_count}"
    else
        printf "FAIL: observed %d PETSc state resets\n" "${reset_count}"
        failures=$((failures + 1))
    fi

    if grep -q "DIVERGED_FUNCTION_DOMAIN" "${solver_log}"; then
        echo "FAIL: stale DIVERGED_FUNCTION_DOMAIN was reported on retry"
        failures=$((failures + 1))
    else
        echo "PASS: no stale DIVERGED_FUNCTION_DOMAIN on retry"
    fi

    if grep -Eq "Segmentation Violation|MPI_ABORT|Caught signal" "${solver_log}"; then
        echo "FAIL: PETSc retry path crashed"
        failures=$((failures + 1))
    else
        echo "PASS: PETSc retry path did not crash"
    fi

    if [[ "${run_status}" -eq 0 ]]; then
        echo "INFO: Allrun returned 0 after the intentional solver abort"
    else
        echo "INFO: Allrun returned ${run_status} after the intentional solver abort"
    fi

    return "${failures}"
}

run_petsc_retry_recovery_check() {
    local run_status=0
    local retry_log="${PETSC_RETRY_RECOVERY_CASE_DIR}/${PETSC_RETRY_RECOVERY_LOGFILE}"
    local solver_log="${PETSC_RETRY_RECOVERY_CASE_DIR}/log.solids4Foam"
    local retry_count=0
    local same_deltaT_retry_count=0
    local reduced_deltaT_retry_count=0
    local max_it_count=0
    local reset_count=0
    local failures=0

    echo "============================================================"
    echo "PETSc retry recovery regression check"
    echo "============================================================"

    if [[ -z "${PETSC_DIR:-}" ]]; then
        echo "Skipping PETSc retry recovery check because PETSc is not installed"
        echo "Please set the PETSC_DIR variable and rebuild solids4foam"
        return 0
    fi

    prepare_case "${PETSC_RETRY_RECOVERY_CASE_DIR}"
    configure_petsc_retry_recovery_case

    ( cd "${PETSC_RETRY_RECOVERY_CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

    set +e
    ( cd "${PETSC_RETRY_RECOVERY_CASE_DIR}" && ./Allrun monolithic > "${PETSC_RETRY_RECOVERY_LOGFILE}" 2>&1 )
    run_status=$?
    set -e

    if declare -F "solids4Foam::regressionCaseSkipped" > /dev/null \
        && solids4Foam::regressionCaseSkipped "${retry_log}"
    then
        echo "Skipping PETSc retry recovery check because the tutorial skipped in this environment"
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
        echo "PASS: retry recovery path reached End"
    else
        echo "FAIL: retry recovery path did not reach End"
        failures=$((failures + 1))
    fi

    if grep -Eq "Segmentation Violation|MPI_ABORT|Caught signal|^ERROR$|FOAM FATAL ERROR|Exceeded the maximum number|minimum allowed time" "${solver_log}"; then
        echo "FAIL: PETSc retry recovery path crashed or aborted"
        failures=$((failures + 1))
    else
        echo "PASS: PETSc retry recovery path did not crash"
    fi

    return "${failures}"
}

# ------------------------------------------------------------
# Clean & run case
# ------------------------------------------------------------

CHECK_ONLY=false
PETSC_RETRY_ONLY=false
PETSC_RETRY_RECOVERY_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run)
            CHECK_ONLY=true
            ;;
        --petsc-retry-only)
            PETSC_RETRY_ONLY=true
            ;;
        --petsc-retry-recovery-only)
            PETSC_RETRY_RECOVERY_ONLY=true
            ;;
        *)
            ;;
    esac
done

if [ "${PETSC_RETRY_ONLY}" = true ]; then
    run_petsc_retry_check
    exit $?
fi

if [ "${PETSC_RETRY_RECOVERY_ONLY}" = true ]; then
    run_petsc_retry_recovery_check
    exit $?
fi

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# OpenFOAM variant compatibility
mkdir -p "${CASE_DIR}/postProcessing/fluid/forces/0"
(
    cd "${CASE_DIR}/postProcessing/fluid/forces/0"

    # foam-extend writes forces to a 'forces' sub-directory so we will create a
    # link
    if [[ ! -e force.dat && -f ../../../../forces/0/forces.dat ]]; then
        ln -s ../../../../forces/0/forces.dat force.dat
    fi

    # OpenFOAM.org creates forces.dat instead of force.dat
    if [[ ! -e force.dat && -f forces.dat ]]; then
        ln -s forces.dat force.dat
    fi
)

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_vertical_displacement() {
    awk '{print $3}' "${CASE_DIR}/${DISP_FILE}" | tail -1
}

extract_mean_force_tail() {
    tail -n "${FORCE_AVG_SAMPLES}" "${CASE_DIR}/${FORCE_FILE}" | \
    awk '
    {
        # Remove parentheses
        gsub(/[()]/, "", $0)

        # After cleanup, fields are:
        # $1 = time
        # $3 = Fy (both formats)
        sum += $3
        n++
    }
    END {
        if (n > 0) print sum/n
    }'
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

max_disp=$(extract_vertical_displacement)
mean_force=$(extract_mean_force_tail)

if [[ -z "${max_disp}" || -z "${mean_force}" ]]; then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

disp_diff=$(awk "BEGIN {print ${max_disp} - ${REF_MAX_DISP}}")
disp_diff_abs=$(abs "${disp_diff}")

force_diff=$(awk "BEGIN {print ${mean_force} - ${REF_MEAN_FORCE}}")
force_diff_abs=$(abs "${force_diff}")

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${disp_diff_abs} < ${DISP_MAX_TOL})}"; then
    printf "PASS: max displacement = %.6g (Δ = %.3g)\n" \
        "${max_disp}" "${disp_diff_abs}"
else
    printf "FAIL: max displacement = %.6g (Δ = %.3g)\n" \
        "${max_disp}" "${disp_diff_abs}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${force_diff_abs} < ${FORCE_MEAN_TOL})}"; then
    printf "PASS: mean force = %.6g (Δ = %.3g)\n" \
        "${mean_force}" "${force_diff_abs}"
else
    printf "FAIL: mean force = %.6g (Δ = %.3g)\n" \
        "${mean_force}" "${force_diff_abs}"
    failures=$((failures + 1))
fi

if [ "$CHECK_ONLY" = false ]; then
    if ! run_petsc_retry_check; then
        failures=$((failures + 1))
    fi

    if ! run_petsc_retry_recovery_check; then
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
