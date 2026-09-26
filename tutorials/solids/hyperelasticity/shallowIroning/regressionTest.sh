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
# shallowIroning regression test
# The full ironing stroke (endTime 6) takes too long for a regression test,
# so only the first half of the vertical indentation stage is simulated, up
# to t = 0.5 s (the block comes into contact with the foundation at
# t = 0.21 s). The checks are:
#   1. the momentum loop converged, i.e. did not stop at the maximum number
#      of correctors, in every time step of the shortened run, so that the
#      checked values do not depend on the iteration path;
#   2. the total reaction force on the "displacement" patch (top of the
#      block) at t = 0.25 s and t = 0.5 s against the values this test was
#      calibrated with (regression check);
#   3. the vertical reaction force at t = 0.5 s lies within the spread of
#      the reference solutions in the reference directory at that instant,
#      linearly interpolated (benchmark check).
# ============================================================

REGRESSION_END_TIME=0.5
EARLY_CHECK_TIME=0.25

# Calibrated with foam-extend-4.1, with the material-aware gradient that two
# materials need on the constitutive law framework (leastSquaresS4f; the
# removed legacy model ran extendedLeastSquares and gave -10.7944, 15.6342 and
# -90.238):
#   force_y = -10.8592 N at t = 0.25 s
#   force_x = 14.8288 N and force_y = -90.0582 N at t = 0.5 s
# The bands are +/- 1% (force_y) and +/- 2% (force_x).
EARLY_FORCE_Y_MIN=-10.97
EARLY_FORCE_Y_MAX=-10.75
FORCE_X_MIN=14.53
FORCE_X_MAX=15.13
FORCE_Y_MIN=-90.96
FORCE_Y_MAX=-89.16

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
FORCE_FILE="postProcessing/0/solidForcesdisplacement.dat"
REFERENCE_VERTICAL_FILES=(
    reference/asterVertical.dat
    reference/fischerWriggersVertical.dat
    reference/hartmannOliverVertical.dat
    reference/pouliosRenardVertical.dat
)

echo "============================================================"
echo "shallowIroning regression test (up to t = ${REGRESSION_END_TIME} s)"
echo "Momentum loop converged in every time step"
echo "force_y at t = ${EARLY_CHECK_TIME} s in [${EARLY_FORCE_Y_MIN}, ${EARLY_FORCE_Y_MAX}] N"
echo "force_x at t = ${REGRESSION_END_TIME} s in [${FORCE_X_MIN}, ${FORCE_X_MAX}] N"
echo "force_y at t = ${REGRESSION_END_TIME} s in [${FORCE_Y_MIN}, ${FORCE_Y_MAX}] N"
echo "Vertical force at t = ${REGRESSION_END_TIME} s within the reference solutions"
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

    sed -i.bak "s/^endTime[[:space:]]\+6;/endTime         ${REGRESSION_END_TIME};/" \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"

    if ! grep -q "^endTime[[:space:]]\+${REGRESSION_END_TIME};" \
        "${CASE_DIR}/system/controlDict"
    then
        echo "FAIL: Could not set endTime to ${REGRESSION_END_TIME} in system/controlDict"
        exit 1
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

failures=0

if [[ ! -f "${CASE_DIR}/${SOLVER_LOGFILE}" ]]; then
    echo "FAIL: Could not find ${SOLVER_LOGFILE}"
    exit 1
fi

if grep -qE 'FOAM FATAL|^ERROR$|\[stack trace\]' "${CASE_DIR}/${SOLVER_LOGFILE}" \
    || ! grep -q '^End' "${CASE_DIR}/${SOLVER_LOGFILE}"
then
    echo "FAIL: solids4Foam did not complete"
    exit 1
fi

if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
    echo "FAIL: Could not find ${FORCE_FILE}"
    exit 1
fi

# Value in the force history file at a given time (to within half a step)
value_at() {
    local column="$1" time="$2"
    awk -v c="${column}" -v t="${time}" \
        '!/^#/ && ($1 - t)^2 < 0.000025 {v = $c} END {print v}' \
        "${CASE_DIR}/${FORCE_FILE}"
}

# Linear interpolation of a reference file at a given (normalised) time,
# scaled by 100 to give the force (N) on the 1 mm thick domain
reference_at() {
    local file="$1" time="$2"
    awk -v x="${time}" \
        '!/^#/ && NF >= 2 {
            if (n && px <= x && $1 >= x) {
                print 100*(py + ($2 - py)*(x - px)/($1 - px)); found = 1; exit
            }
            px = $1; py = $2; n = 1
        }
        END {if (!found) exit 1}' "${file}"
}

# Check that the momentum loop converged in the time step ending at a given
# time, i.e. it did not stop at the maximum number of correctors
check_converged() {
    local time="$1"
    local status
    status=$(awk -v t="${time}" '
        /^Time = / {cur = $3; next}
        (cur - t)^2 < 1e-8 && /Max iterations reached/ {s = "capped"}
        (cur - t)^2 < 1e-8 && /residual has converged|Both residuals have converged/ {s = "converged"}
        END {print s}' "${CASE_DIR}/${SOLVER_LOGFILE}")

    if [[ "${status}" == "converged" ]]; then
        printf "PASS: Momentum loop converged at t = %s\n" "${time}"
    else
        printf "FAIL: Momentum loop did not converge at t = %s\n" "${time}"
        failures=$((failures + 1))
    fi
}

# Check that no time step reached the maximum number of correctors
check_no_capped_steps() {
    local nCapped
    nCapped=$(grep -c "Max iterations reached" "${CASE_DIR}/${SOLVER_LOGFILE}" || true)

    if [[ "${nCapped}" == "0" ]]; then
        echo "PASS: No time step reached the maximum number of correctors"
    else
        echo "FAIL: ${nCapped} time step(s) reached the maximum number of correctors"
        failures=$((failures + 1))
    fi
}

# Check that a value lies within [min, max]
check_range() {
    local label="$1" value="$2" min="$3" max="$4"
    if awk "BEGIN {exit !(${value} >= ${min} && ${value} <= ${max})}"; then
        printf "PASS: %s = %.6g, expected [%s, %s]\n" "${label}" "${value}" "${min}" "${max}"
    else
        printf "FAIL: %s = %.6g, expected [%s, %s]\n" "${label}" "${value}" "${min}" "${max}"
        failures=$((failures + 1))
    fi
}

early_force_y=$(value_at 3 "${EARLY_CHECK_TIME}")
force_x=$(value_at 2 "${REGRESSION_END_TIME}")
force_y=$(value_at 3 "${REGRESSION_END_TIME}")

if [[ -z "${early_force_y}" || -z "${force_x}" || -z "${force_y}" ]]; then
    echo "FAIL: Could not extract the force at t = ${EARLY_CHECK_TIME} and t = ${REGRESSION_END_TIME}"
    exit 1
fi

check_no_capped_steps
check_converged "${EARLY_CHECK_TIME}"
check_converged "${REGRESSION_END_TIME}"

check_range "force_y at t = ${EARLY_CHECK_TIME}" "${early_force_y}" \
    "${EARLY_FORCE_Y_MIN}" "${EARLY_FORCE_Y_MAX}"
check_range "force_x at t = ${REGRESSION_END_TIME}" "${force_x}" \
    "${FORCE_X_MIN}" "${FORCE_X_MAX}"
check_range "force_y at t = ${REGRESSION_END_TIME}" "${force_y}" \
    "${FORCE_Y_MIN}" "${FORCE_Y_MAX}"

# Spread of the reference vertical forces at the check time; the
# normalised time equals the solids4foam time during the indentation stage
ref_min=""
ref_max=""
for file in "${REFERENCE_VERTICAL_FILES[@]}"; do
    ref=$(reference_at "${SCRIPT_DIR}/${file}" "${REGRESSION_END_TIME}")
    ref_min=$(awk -v a="${ref_min}" -v b="${ref}" 'BEGIN {print (a == "" || b < a) ? b : a}')
    ref_max=$(awk -v a="${ref_max}" -v b="${ref}" 'BEGIN {print (a == "" || b > a) ? b : a}')
done

# The force on the displacement patch points upwards (negative y)
check_range "Vertical force at t = ${REGRESSION_END_TIME} vs reference solutions" \
    "$(awk "BEGIN {print -(${force_y})}")" "${ref_min}" "${ref_max}"

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
