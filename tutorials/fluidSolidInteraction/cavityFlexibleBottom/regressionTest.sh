#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

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
# Reference values updated for the interface-normal correction (PR #375): the
# fluid pressure is now applied using the deformed interface normals rather than
# the initial-configuration ones, which shifts every FSI result. See
# https://github.com/solids4foam/solids4foam/pull/375
# Monitor point moved from (4 -1 0.5), which snapped to the plate underside, to
# point A on the top surface, (4 0 0.5); reference regenerated accordingly.
REF_MAX_DISP=-0.2045
REF_MEAN_FORCE=-5.150

# The quasi-monolithic Newton coupling, which needs OpenFOAM.com and PETSc, is
# checked against the same references: its damped solid settles to the same
# steady state, to 7e-4 in displacement and 0.03 in mean force

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

echo "============================================================"
echo "cavityFlexibleBottom FSI regression test"
echo "Max displacement difference < ${DISP_MAX_TOL}"
echo "Mean force difference       < ${FORCE_MEAN_TOL}"
echo "============================================================"
echo

# The helpers below act on CASE_DIR, which run_and_check sets for each arm

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

    # Copy the record of a selected case variant too, so that the Allrun of
    # the copy can restore the default variant
    if [[ -f "${SCRIPT_DIR}/.solids4FoamCaseVariant" ]]; then
        cp -a "${SCRIPT_DIR}/.solids4FoamCaseVariant" "${CASE_DIR}/"
    fi
}

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_vertical_displacement() {
    awk '{print $3}' "${CASE_DIR}/${DISP_FILE}" | tail -1
}

extract_mean_force_tail() {
    # The forces functionObject writes a different set of columns depending on
    # the OpenFOAM version: OpenFOAM.com writes the total force followed by the
    # pressure and viscous contributions, whereas OpenFOAM.org and foam-extend
    # write the pressure and viscous contributions followed by the moments. The
    # number of columns is used to tell them apart, so that the total force is
    # compared in both cases.
    tail -n "${FORCE_AVG_SAMPLES}" "${CASE_DIR}/${FORCE_FILE}" | \
    awk '
    {
        # Remove parentheses
        gsub(/[()]/, "", $0)

        if (NF >= 13)
        {
            # time, pressure, viscous, moments: sum the contributions
            sum += $3 + $6
        }
        else
        {
            # time, total, pressure, viscous: use the total directly
            sum += $3
        }
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
# Run and check one arm
# Arguments: arm name and Allrun arguments
# ------------------------------------------------------------

run_and_check() {
    local ARM="$1"
    local ALLRUN_ARGS="$2"
    local CASE_DIR="${REGRESSION_ROOT}/${ARM}"

    echo "Running the ${ARM} arm"

    if [ "$CHECK_ONLY" = false ]; then
        prepare_case
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
        ( cd "${CASE_DIR}" && ./Allrun ${ALLRUN_ARGS} > "${ALLRUN_LOGFILE}" 2>&1 )
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
    # Extract values
    # ------------------------------------------------------------

    max_disp=$(extract_vertical_displacement)
    mean_force=$(extract_mean_force_tail)

    if [[ -z "${max_disp}" || -z "${mean_force}" ]]; then
        echo "FAIL [${ARM}]: Could not extract regression quantities"
        failures=$((failures + 1))
        return
    fi

    disp_diff=$(awk "BEGIN {print ${max_disp} - ${REF_MAX_DISP}}")
    disp_diff_abs=$(abs "${disp_diff}")

    force_diff=$(awk "BEGIN {print ${mean_force} - ${REF_MEAN_FORCE}}")
    force_diff_abs=$(abs "${force_diff}")

    # ------------------------------------------------------------
    # Checks
    # ------------------------------------------------------------

    if awk "BEGIN {exit !(${disp_diff_abs} < ${DISP_MAX_TOL})}"; then
        printf "PASS [%s]: max displacement = %.6g (Δ = %.3g)\n" \
            "${ARM}" "${max_disp}" "${disp_diff_abs}"
    else
        printf "FAIL [%s]: max displacement = %.6g (Δ = %.3g)\n" \
            "${ARM}" "${max_disp}" "${disp_diff_abs}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${force_diff_abs} < ${FORCE_MEAN_TOL})}"; then
        printf "PASS [%s]: mean force = %.6g (Δ = %.3g)\n" \
            "${ARM}" "${mean_force}" "${force_diff_abs}"
    else
        printf "FAIL [%s]: mean force = %.6g (Δ = %.3g)\n" \
            "${ARM}" "${mean_force}" "${force_diff_abs}"
        failures=$((failures + 1))
    fi
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

failures=0

run_and_check partitioned ""

# The quasi-monolithic coupling needs OpenFOAM.com and PETSc
if [[ "${WM_PROJECT_VERSION:-}" == v* && -n "${PETSC_DIR:-}" ]]; then
    run_and_check monolithic monolithic
else
    echo "Skipping the monolithic arm: it needs OpenFOAM.com and PETSc"
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
