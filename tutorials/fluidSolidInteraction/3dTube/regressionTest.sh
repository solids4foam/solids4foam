#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
BACKWARD_CASE_DIR="${REGRESSION_ROOT}/backwardRestart"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# Beam-in-cross-flow FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

DISP_MAX_TOL=1e-6       # max displacement absolute tolerance
FORCE_MEAN_TOL=2.6e-2   # mean force tolerance

# Regression end time for the copied case only
REG_END_TIME=0.0015
BACKWARD_END_TIME=5e-5
BACKWARD_WRITE_INTERVAL=2.5e-5

# Number of samples from end of force.dat to average
FORCE_AVG_SAMPLES=50

# Reference values at REG_END_TIME
REF_MAX_DISP=2.23646e-07
REF_MEAN_FORCE=0.0320942

# The final probe displacement magnitude of the removed legacy
# mechanicalModel, from the last commit that had it (mcl-stage8-coverage,
# c3a92b3d), per fork: the coupled answer moves between forks by more than the
# framework moved from legacy on any one of them. The material is linear
# elastic, so the framework solves the same problem, and it reproduced these in
# all eight digits printed. The comparison it replaces allowed 1e-6 absolute,
# which is more than the value itself and so could not fail; 1e-6 relative is
# the same agreement stated in a form that can
case "$(solids4Foam::foamFlavour)" in
    com)        LEGACY_FINAL_MAGD=7.18043e-07 ;;
    org)        LEGACY_FINAL_MAGD=7.06215e-07 ;;
    foamextend) LEGACY_FINAL_MAGD=7.08637e-07 ;;
esac
LEGACY_MAGD_REL_TOL=1e-6

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_displacement.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

echo "============================================================"
echo "3dTube FSI regression test"
echo "Regression end time         = ${REG_END_TIME}"
echo "Max displacement difference < ${DISP_MAX_TOL}"
echo "Mean force difference       < ${FORCE_MEAN_TOL}"
echo "============================================================"
echo

copy_case() {
    local destination="$1"

    rm -rf "${destination}"
    mkdir -p "${destination}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${destination}/"
    done
}

prepare_case() {
    copy_case "${CASE_DIR}"

    sed -i "s/^\(endTime[[:space:]]*\).*/\1${REG_END_TIME};/" "${CASE_DIR}/system/controlDict"
}

# Against the legacy answer.
#
# The first framework coverage of a fluid-solid interaction case. It matters
# beyond this tutorial: elasticWallPressure looks impK up from the solid mesh
# by name, to build the p-wave speed it uses for its added-mass term, so
# frameworkImpK() has to register a field the fluid side can find. The same
# lookup is used by the contact penalty models and the cohesive zone models
check_against_legacy() {
    local failures=0

    local solver_log
    solver_log=$(find "${CASE_DIR}" -name 'log.solids4Foam' | tail -n 1)

    if [[ -n "${solver_log}" ]] \
        && grep -q "Selecting mechanical constitutive law" "${solver_log}"
    then
        echo "PASS: the solid took its material from the framework"
    else
        echo "FAIL: the solver log shows no mechanical constitutive law"
        failures=$((failures + 1))
    fi

    # Column 5 is magD. Column 2 is Dx, which is identically zero on this
    # geometry, so comparing it would compare zero with zero and pass whatever
    # the solver did
    local disp
    disp=$(awk 'END {print $5}' "${CASE_DIR}/${DISP_FILE}")

    # A comparison against zero is vacuous: it would pass for a run that
    # produced nothing at all
    if [[ -z "${disp}" ]] || awk "BEGIN {exit !(${disp} <= 0)}"; then
        echo "FAIL: the final displacement is '${disp}', so this comparison"
        echo "      would pass whatever the solver did"
        return $((failures + 1))
    fi

    if awk "BEGIN {d = ${disp} - ${LEGACY_FINAL_MAGD}; if (d < 0) d = -d;
                   exit !(d <= ${LEGACY_MAGD_REL_TOL} * ${LEGACY_FINAL_MAGD})}"
    then
        printf "PASS: final displacement matches the legacy model (%.8g vs %.8g)\n" \
            "${disp}" "${LEGACY_FINAL_MAGD}"
    else
        printf "FAIL: final displacement differs from the legacy model (%.8g vs %.8g)\n" \
            "${disp}" "${LEGACY_FINAL_MAGD}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

prepare_backward_case() {
    copy_case "${BACKWARD_CASE_DIR}"

    # This total-strain material does not need constitutive kinematic history,
    # but a restart must say so explicitly.
    sed -i \
        's/^    nCorrectors/    restart                 no;\n\n    nCorrectors/' \
        "${BACKWARD_CASE_DIR}/constant/solid/solidProperties"
    sed -i "s/^\(endTime[[:space:]]*\).*/\1${BACKWARD_END_TIME};/" \
        "${BACKWARD_CASE_DIR}/system/controlDict"
    sed -i "s/^\(writeInterval[[:space:]]*\).*/\1${BACKWARD_WRITE_INTERVAL};/" \
        "${BACKWARD_CASE_DIR}/system/controlDict"
    sed -i "s/^\(startFrom[[:space:]]*\).*/\1latestTime;/" \
        "${BACKWARD_CASE_DIR}/system/controlDict"
    sed -i "s/default[[:space:]]*Euler;/default            backward;/" \
        "${BACKWARD_CASE_DIR}/system/fluid/fvSchemes"
}

run_backward_restart_test() {
    prepare_backward_case
    (
        cd "${BACKWARD_CASE_DIR}"
        ./Allclean > /dev/null 2>&1 || true
        ./Allrun > log.Allrun 2>&1
    )
}

check_backward_restart() {
    (
        cd "${BACKWARD_CASE_DIR}"
        Test-fluxCorrectedVelocityRestart \
            > log.Test-fluxCorrectedVelocityRestart 2>&1
    )
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

find_force_file() {
    local candidate
    for candidate in \
        "${CASE_DIR}/postProcessing/fluid/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/fluid/forces/0/forces.dat" \
        "${CASE_DIR}/postProcessing/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/forces/0/forces.dat" \
        "${CASE_DIR}/postProcessing/0/solidForcesinner-wall.dat" \
        "${CASE_DIR}/postProcessing/0/solidForcesDisplacementsloading.dat"
    do
        if [[ -f "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done
    return 1
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

disp_time=$(latest_numeric_time "${CASE_DIR}/${DISP_FILE}" || true)
force_file=""
if force_file=$(find_force_file); then
    force_time=$(latest_numeric_time "${force_file}" || true)
else
    force_time=""
fi

# The one supported reason not to run is the environment the case declares it
# needs (PETSc); Allrun then says so, and that is a skip. A run that starts
# and stops short is a failure, not a skip: it used to exit 0 here, which kept
# a crash on OpenFOAM.org and foam-extend green for as long as it lasted (#455)
if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the case does not run in this environment"
    exit 0
fi

if [[ -z "${disp_time}" || -z "${force_file}" || -z "${force_time}" ]]; then
    echo "FAIL: the case produced no displacement or force history"
    exit 1
fi

if ! awk "BEGIN {exit !(${disp_time} + 0 >= ${REG_END_TIME})}"; then
    echo "FAIL: the case stopped at ${disp_time}, before the requested end time ${REG_END_TIME}"
    exit 1
fi

if ! awk "BEGIN {exit !(${force_time} + 0 >= ${REG_END_TIME})}"; then
    echo "FAIL: the force history stopped at ${force_time}, before the requested end time ${REG_END_TIME}"
    exit 1
fi

if [ "$CHECK_ONLY" = false ]; then
    run_backward_restart_test
fi
check_backward_restart

# Said, so that a run which reached this check can be told from one which
# skipped it: the check used to be skipped silently on two forks (#455).
# A failure has already stopped the script under set -e
echo "PASS: fluxCorrectedVelocity backward restart" \
    "($(grep -o 'Checked [0-9]* fluxCorrectedVelocity patches' \
        "${BACKWARD_CASE_DIR}/log.Test-fluxCorrectedVelocityRestart"))"

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

extract_max_displacement() {
    awk 'NR > 1 {print $3}' "${CASE_DIR}/${DISP_FILE}" | sort -g | tail -1
}

extract_mean_force_tail() {
    tail -n "${FORCE_AVG_SAMPLES}" "$1" | \
    awk '
    {
        gsub(/[()]/, "", $0)
        sum += $2
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

max_disp=$(extract_max_displacement)
mean_force=$(extract_mean_force_tail "${force_file}")

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

check_against_legacy || failures=$((failures + $?))

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
