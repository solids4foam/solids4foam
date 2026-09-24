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
# twistingHemisphere regression test
# To keep the run short, only the start of the load history is simulated:
# the full indentation (t = 0 to 4) and the first 40 degrees of twist
# (t = 4 to 6; angle = 20*(t - 4) degrees).
#
# The checks are made only at time steps whose momentum loop converged.
# On the shipped coarse mesh, the increments from t = 3.8 to 5.7 (the end of
# the indentation and the stick-to-slip transition at the start of the
# twist) reach the 2000-corrector limit without meeting the tolerance, so
# their values depend on the iteration path and are not used. The checks are:
#   1. the vertical force at t = 3.5 (indentation) and the vertical force
#      and twisting torque at t = 6 (40 degrees) against the values this
#      test was calibrated with (regression check);
#   2. the vertical force and twisting torque at t = 6 against the digitised
#      solution of Sauer and De Lorenzis at 40 degrees, linearly
#      interpolated from reference/deLorenzisForce.dat and
#      reference/deLorenzisMoment.dat (benchmark check). The torque tolerance is wide because, on the
#      coarse mesh, the torque rises more slowly than in the reference.
# ============================================================

REGRESSION_END_TIME=6
REGRESSION_ANGLE=40
INDENTATION_CHECK_TIME=3.5

# Calibrated with foam-extend-4.1 on the shipped coarse mesh:
#   vertical force = 1.01883 at t = 3.5
#   vertical force = 1.23495 and twisting torque = 0.218597 at t = 6
INDENTATION_FORCE_MIN=1.008
INDENTATION_FORCE_MAX=1.029
FORCE_MIN=1.225
FORCE_MAX=1.245
TORQUE_MIN=0.212
TORQUE_MAX=0.225

# Allowed relative deviations from Sauer and De Lorenzis at 40 degrees
# (calibration run: force -0.8%, torque -18.0%)
REFERENCE_FORCE_TOLERANCE=0.03
REFERENCE_TORQUE_TOLERANCE=0.25

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
FORCE_FILE="postProcessing/0/solidForcessphere-displacement.dat"
TORQUE_FILE="postProcessing/0/solidTorquesphere-displacementsphereTorque.dat"
REFERENCE_FORCE_FILE="reference/deLorenzisForce.dat"
REFERENCE_TORQUE_FILE="reference/deLorenzisMoment.dat"

echo "============================================================"
echo "twistingHemisphere regression test (indentation + ${REGRESSION_ANGLE} degrees of twist)"
echo "Vertical force at t = ${INDENTATION_CHECK_TIME} in [${INDENTATION_FORCE_MIN}, ${INDENTATION_FORCE_MAX}]"
echo "Vertical force at t = ${REGRESSION_END_TIME} in [${FORCE_MIN}, ${FORCE_MAX}]"
echo "Twisting torque at t = ${REGRESSION_END_TIME} in [${TORQUE_MIN}, ${TORQUE_MAX}]"
echo "Force and torque at ${REGRESSION_ANGLE} degrees within" \
     "$(awk "BEGIN {print 100*${REFERENCE_FORCE_TOLERANCE}}")% and" \
     "$(awk "BEGIN {print 100*${REFERENCE_TORQUE_TOLERANCE}}")% of Sauer and De Lorenzis"
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

    sed -i.bak "s/^endTime[[:space:]]\+13;/endTime         ${REGRESSION_END_TIME};/" \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"
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

if grep -qE "FOAM FATAL|^ERROR$|stack trace" "${CASE_DIR}/${SOLVER_LOGFILE}"; then
    echo "FAIL: ${SOLVER_LOGFILE} reports an error"
    failures=$((failures + 1))
fi

for file in "${FORCE_FILE}" "${TORQUE_FILE}"; do
    if [[ ! -f "${CASE_DIR}/${file}" ]]; then
        echo "FAIL: Could not find ${file}"
        exit 1
    fi
done

# Value in a history file at a given time (to within half a time step)
value_at() {
    local file="$1" column="$2" time="$3"
    awk -v c="${column}" -v t="${time}" \
        '!/^#/ && ($1 - t)^2 < 0.0025 {v = $c} END {print v}' "${file}"
}

# Linear interpolation of a comma-separated reference file at a given angle
reference_at() {
    local file="$1" angle="$2"
    awk -F'[, \t]+' -v x="${angle}" \
        '!/^#/ && NF >= 2 {
            if (n && px <= x && $1 >= x) {
                print py + ($2 - py)*(x - px)/($1 - px); found = 1; exit
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

# Check that a value lies within [min, max]
check_range() {
    local label="$1" value="$2" min="$3" max="$4"
    if awk "BEGIN {exit !(${value} >= ${min} && ${value} <= ${max})}"; then
        printf "PASS: %s: %.6g\n" "${label}" "${value}"
    else
        printf "FAIL: %s: %.6g\n" "${label}" "${value}"
        failures=$((failures + 1))
    fi
}

# Check the relative deviation from the reference solution
check_reference() {
    local label="$1" value="$2" reference="$3" tolerance="$4"
    local deviation
    deviation=$(awk "BEGIN {d = ${value}/${reference} - 1; print (d < 0 ? -d : d)}")
    if awk "BEGIN {exit !(${deviation} <= ${tolerance})}"; then
        printf "PASS: %s: %.6g vs %.6g (%.2f%% from Sauer and De Lorenzis)\n" \
            "${label}" "${value}" "${reference}" "$(awk "BEGIN {print 100*${deviation}}")"
    else
        printf "FAIL: %s: %.6g vs %.6g (%.2f%% from Sauer and De Lorenzis)\n" \
            "${label}" "${value}" "${reference}" "$(awk "BEGIN {print 100*${deviation}}")"
        failures=$((failures + 1))
    fi
}

indentation_force=$(value_at "${CASE_DIR}/${FORCE_FILE}" 3 "${INDENTATION_CHECK_TIME}")
force=$(value_at "${CASE_DIR}/${FORCE_FILE}" 3 "${REGRESSION_END_TIME}")
torque=$(value_at "${CASE_DIR}/${TORQUE_FILE}" 2 "${REGRESSION_END_TIME}")

if [[ -z "${indentation_force}" || -z "${force}" || -z "${torque}" ]]; then
    echo "FAIL: Could not extract the force and torque at t = ${INDENTATION_CHECK_TIME} and t = ${REGRESSION_END_TIME}"
    exit 1
fi

reference_force=$(reference_at "${SCRIPT_DIR}/${REFERENCE_FORCE_FILE}" "${REGRESSION_ANGLE}")
reference_torque=$(reference_at "${SCRIPT_DIR}/${REFERENCE_TORQUE_FILE}" "${REGRESSION_ANGLE}")

# The force on the sphere-displacement patch points upwards (negative y)
indentation_force=$(awk "BEGIN {print -(${indentation_force})}")
force=$(awk "BEGIN {print -(${force})}")

check_converged "${INDENTATION_CHECK_TIME}"
check_converged "${REGRESSION_END_TIME}"

check_range "Vertical force at t = ${INDENTATION_CHECK_TIME}" "${indentation_force}" \
    "${INDENTATION_FORCE_MIN}" "${INDENTATION_FORCE_MAX}"
check_range "Vertical force at t = ${REGRESSION_END_TIME}" "${force}" \
    "${FORCE_MIN}" "${FORCE_MAX}"
check_range "Twisting torque at t = ${REGRESSION_END_TIME}" "${torque}" \
    "${TORQUE_MIN}" "${TORQUE_MAX}"

check_reference "Vertical force at ${REGRESSION_ANGLE} degrees" "${force}" \
    "${reference_force}" "${REFERENCE_FORCE_TOLERANCE}"
check_reference "Twisting torque at ${REGRESSION_ANGLE} degrees" "${torque}" \
    "${reference_torque}" "${REFERENCE_TORQUE_TOLERANCE}"

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
