#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# blobInTreacle FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression tolerances (relative)
# ------------------------------------------------------------

# The steady displacement and force converge tightly; the displacement at
# t = 1 s is taken during the ramp, where the partitioned coupling and the
# added-mass effect matter most
DISP_T1_TOL=0.01
DISP_END_TOL=0.005
FORCE_END_TOL=0.01

# Reference values: OpenFOAM v2412, serial, Apple M1 Ultra
REF_DISP_T1=0.0565316
REF_DISP_END=0.117228
REF_FORCE_END=15.8217

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

echo "============================================================"
echo "blobInTreacle FSI regression test"
echo "x displacement at t = 1 s relative difference < ${DISP_T1_TOL}"
echo "Final x displacement relative difference      < ${DISP_END_TOL}"
echo "Final x force relative difference             < ${FORCE_END_TOL}"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" \
           || "${base_item}" == "verification" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done
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

extract_displacement_at() {
    awk -v t="$1" '
        $1 !~ /^#/ && ($1 - t < 1e-8 && t - $1 < 1e-8) { print $2 }
    ' "${CASE_DIR}/${DISP_FILE}" | tail -1
}

extract_final_displacement() {
    awk '$1 !~ /^#/ { value = $2 } END { print value }' \
        "${CASE_DIR}/${DISP_FILE}"
}

extract_final_force() {
    # The forces functionObject writes a different set of columns depending on
    # the OpenFOAM version: OpenFOAM.com writes the total force followed by the
    # pressure and viscous contributions, whereas OpenFOAM.org and foam-extend
    # write the pressure and viscous contributions followed by the moments. The
    # number of columns is used to tell them apart, so that the total force is
    # compared in both cases.
    tail -n 1 "${CASE_DIR}/${FORCE_FILE}" | \
    awk '
    {
        # Remove parentheses
        gsub(/[()]/, "", $0)

        if (NF >= 13)
        {
            # time, pressure, viscous, moments: sum the contributions
            print $2 + $5
        }
        else
        {
            # time, total, pressure, viscous: use the total directly
            print $2
        }
    }'
}

relative_difference() {
    awk -v a="$1" -v b="$2" 'BEGIN {
        d = (a - b)/b
        print (d < 0 ? -d : d)
    }'
}

# ------------------------------------------------------------
# Check that the solver completed
# ------------------------------------------------------------

# Allrun does not return the solver status, so check the solver log and the
# time of the last displacement sample
SOLVER_LOG="${CASE_DIR}/log.solids4Foam"
if [[ ! -f "${SOLVER_LOG}" ]] || ! grep -q "^End" "${SOLVER_LOG}" \
    || grep -q "FOAM FATAL" "${SOLVER_LOG}"; then
    echo "FAIL: solids4Foam did not run to completion; see ${SOLVER_LOG}"
    exit 1
fi

END_TIME=$(awk '/^endTime/ { gsub(";", "", $2); print $2 }' \
    "${CASE_DIR}/system/controlDict")
last_time=$(awk '$1 !~ /^#/ { value = $1 } END { print value }' \
    "${CASE_DIR}/${DISP_FILE}")
if ! awk -v a="${last_time}" -v b="${END_TIME}" \
    'BEGIN { d = a - b; exit !(d < 1e-8 && d > -1e-8) }'; then
    echo "FAIL: last displacement sample at t = ${last_time}," \
        "not at the end time ${END_TIME}"
    exit 1
fi

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

disp_t1=$(extract_displacement_at 1)
disp_end=$(extract_final_displacement)
force_end=$(extract_final_force)

if [[ -z "${disp_t1}" || -z "${disp_end}" || -z "${force_end}" ]]; then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check() {
    local label="$1"
    local value="$2"
    local reference="$3"
    local tolerance="$4"
    local difference
    difference=$(relative_difference "${value}" "${reference}")

    if awk "BEGIN {exit !(${difference} < ${tolerance})}"; then
        printf "PASS: %s = %.6g (relative difference %.3g)\n" \
            "${label}" "${value}" "${difference}"
    else
        printf "FAIL: %s = %.6g (relative difference %.3g)\n" \
            "${label}" "${value}" "${difference}"
        failures=$((failures + 1))
    fi
}

check "x displacement at t = 1 s" "${disp_t1}" "${REF_DISP_T1}" "${DISP_T1_TOL}"
check "final x displacement" "${disp_end}" "${REF_DISP_END}" "${DISP_END_TOL}"
check "final x force" "${force_end}" "${REF_FORCE_END}" "${FORCE_END_TOL}"

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
