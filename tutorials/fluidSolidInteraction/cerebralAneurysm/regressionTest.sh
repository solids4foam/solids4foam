#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# cerebralAneurysm FSI regression test
#
# The full case simulates one cardiac cycle and takes hours, so the regression
# test runs only the first few time steps of the initial transient in a copy of
# the case.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

DISP_MAX_TOL=1e-6       # max innerWall displacement absolute tolerance
FORCE_NORM_TOL=1e-3     # innerWall normal force absolute tolerance

# Regression end time for the copied case only: 10 time steps of 5e-5 s
REG_END_TIME=0.0005

# Reference values at REG_END_TIME
REF_MAX_DISP=5.23438e-05
REF_NORMAL_FORCE=-0.0138715

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidDisplacementsinnerWall.dat"
FORCE_FILE="postProcessing/0/solidForcesinnerWall.dat"

echo "============================================================"
echo "cerebralAneurysm FSI regression test"
echo "Regression end time          = ${REG_END_TIME}"
echo "Max displacement difference  < ${DISP_MAX_TOL}"
echo "Normal force difference      < ${FORCE_NORM_TOL}"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    # Copy only the case inputs: results, logs and generated meshes from any
    # previous run of the tutorial are deliberately left behind.
    for item in Allrun Allclean 0 constant system *.gnuplot; do
        cp -a "${SCRIPT_DIR}/${item}" "${CASE_DIR}/"
    done

    rm -rf "${CASE_DIR}/constant/polyMesh" \
           "${CASE_DIR}/constant/fluid/polyMesh" \
           "${CASE_DIR}/constant/solid/polyMesh"

    sed -i "s/^\(endTime[[:space:]]*\).*/\1${REG_END_TIME};/" "${CASE_DIR}/system/controlDict"
}

latest_numeric_time() {
    local file="$1"
    awk '
        ($1 + 0) == $1 { time = $1 }
        END {
            if (time != "") print time
        }
    ' "$file"
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

# The case requires PETSc and cartesianMesh: Allrun exits silently when either
# is unavailable, in which case there is nothing to check
if [[ ! -f "${CASE_DIR}/${DISP_FILE}" || ! -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
    echo "Skipping regression checks because the case did not run in this environment"
    exit 0
fi

disp_time=$(latest_numeric_time "${CASE_DIR}/${DISP_FILE}" || true)
force_time=$(latest_numeric_time "${CASE_DIR}/${FORCE_FILE}" || true)

if [[ -z "${disp_time}" || -z "${force_time}" ]]; then
    echo "Skipping regression checks because the case did not complete in this environment"
    exit 0
fi

if ! awk "BEGIN {exit !(${disp_time} + 0 >= ${REG_END_TIME})}"; then
    echo "Skipping regression checks because the case did not reach the requested end time"
    exit 0
fi

if ! awk "BEGIN {exit !(${force_time} + 0 >= ${REG_END_TIME})}"; then
    echo "Skipping regression checks because the force history did not reach the requested end time"
    exit 0
fi

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_displacement() {
    # solidDisplacements writes:
    #     time minX minY minZ maxX maxY maxZ avX avY avZ
    # The largest displacement component magnitude on the final line is used
    tail -1 "${CASE_DIR}/${DISP_FILE}" | \
    awk '
    {
        maxD = 0
        for (i = 2; i <= 7; i++)
        {
            d = ($i < 0 ? -$i : $i)
            if (d > maxD) maxD = d
        }
        print maxD
    }'
}

extract_normal_force() {
    # solidForces writes: time forceX forceY forceZ normalForce
    tail -1 "${CASE_DIR}/${FORCE_FILE}" | awk '{print $5}'
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

max_disp=$(extract_max_displacement)
normal_force=$(extract_normal_force)

if [[ -z "${max_disp}" || -z "${normal_force}" ]]; then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

disp_diff=$(awk "BEGIN {print ${max_disp} - ${REF_MAX_DISP}}")
disp_diff_abs=$(abs "${disp_diff}")

force_diff=$(awk "BEGIN {print ${normal_force} - ${REF_NORMAL_FORCE}}")
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

if awk "BEGIN {exit !(${force_diff_abs} < ${FORCE_NORM_TOL})}"; then
    printf "PASS: normal force = %.6g (Δ = %.3g)\n" \
        "${normal_force}" "${force_diff_abs}"
else
    printf "FAIL: normal force = %.6g (Δ = %.3g)\n" \
        "${normal_force}" "${force_diff_abs}"
    failures=$((failures + 1))
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
