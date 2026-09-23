#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# Hron-Turek FSI3 regression test
# ============================================================

# Shortened regression horizon: keep the case quick while still
# exercising the coupled fluid-solid response after coupling starts.
REG_END_TIME=2.5

# Regression tolerances
DISP_TOL=2e-5
FX_TOL=1e-4
FY_TOL=1e-3

# Reference values at REG_END_TIME
# Note: the force references are the total force. They were previously the
# total force plus the pressure force, as the extraction summed the OpenFOAM.com
# total and pressure columns.
# Reference values updated for the interface-normal correction (PR #375): the
# fluid pressure is now applied using the deformed interface normals rather than
# the initial-configuration ones, which shifts every FSI result. See
# https://github.com/solids4foam/solids4foam/pull/375
# The values are the midpoint of OpenFOAM-v2412, OpenFOAM-v2512 and
# OpenFOAM-9; the spread across those is well inside the tolerances above.
REF_TIP_UY=-0.00031014
REF_FX=-0.0393827
REF_FY=-0.0461165

# foam-extend uses GGI rather than AMI for the interface interpolation and has
# a distinct, repeatable tip displacement and force at the regression end time.
if [[ "${WM_PROJECT:-}" == "foam" ]]; then
    REF_TIP_UY=-0.000169319
    REF_FX=-0.0384236
    REF_FY=-0.0442995
fi

ALLRUN_LOGFILE="log.Allrun"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"

echo "============================================================"
echo "Hron-Turek FSI3 regression test"
echo "Regression end time         = ${REG_END_TIME}"
echo "Tip Uy tolerance            < ${DISP_TOL}"
echo "Final Fx tolerance          < ${FX_TOL}"
echo "Final Fy tolerance          < ${FY_TOL}"
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

    sed -i "s/^\(endTime[[:space:]]*\).*/\1${REG_END_TIME};/" "${CASE_DIR}/system/controlDict"
}

run_case() {
    (
        cd "${CASE_DIR}"
        ./Allclean > /dev/null 2>&1 || true
        ./Allrun > "${ALLRUN_LOGFILE}" 2>&1
    )
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
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

extract_final_tip_uy() {
    awk '
    ($1 + 0) == $1 {
        uy = $3
    }
    END {
        if (uy != "") print uy
    }' "${CASE_DIR}/${DISP_FILE}"
}

extract_final_force_components() {
    # The forces functionObject writes a different set of columns depending on
    # the OpenFOAM version: OpenFOAM.com writes the total force followed by the
    # pressure and viscous contributions, whereas OpenFOAM.org and foam-extend
    # write the pressure and viscous contributions followed by the moments. The
    # number of columns is used to tell them apart.
    awk '
    ($1 + 0) == $1 {
        gsub(/[()]/, "", $0)
        if (NF >= 13)
        {
            # time, pressure, viscous, moments: sum the contributions
            fx = $2 + $5
            fy = $3 + $6
        }
        else
        {
            # time, total, pressure, viscous: use the total directly
            fx = $2
            fy = $3
        }
    }
    END {
        if (fx != "" && fy != "") print fx, fy
    }' "$1"
}

find_force_file() {
    # OpenFOAM.com and OpenFOAM.org write the forces under postProcessing,
    # whereas foam-extend writes them to <case>/forces/<startTime>
    local candidate
    for candidate in \
        "${CASE_DIR}/postProcessing/fluid/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/fluid/forces/0/forces.dat" \
        "${CASE_DIR}/postProcessing/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/forces/0/forces.dat" \
        "${CASE_DIR}/forces/0/forces.dat"
    do
        if [[ -f "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done
    return 1
}

prepare_case
run_case

# A skip is only valid if the tutorial declared one in the Allrun log. Anything
# else that leaves the expected output missing or incomplete is a failure.
if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

tip_time=$(latest_numeric_time "${CASE_DIR}/${DISP_FILE}" || true)

if [[ -z "${tip_time}" ]]; then
    echo "FAIL: the case did not run or did not complete in this environment:"
    echo "      displacement output is missing and the tutorial did not declare a skip"
    echo "      (see ${CASE_DIR}/${ALLRUN_LOGFILE})"
    exit 1
fi

if ! awk "BEGIN {exit !(${tip_time} + 0 >= ${REG_END_TIME})}"; then
    echo "FAIL: the tip displacement history stops at t = ${tip_time}, short of the"
    echo "      requested end time ${REG_END_TIME}: the case did not complete"
    exit 1
fi

if ! force_file=$(find_force_file); then
    echo "FAIL: the case did not run or did not complete in this environment:"
    echo "      force output is missing and the tutorial did not declare a skip"
    echo "      (see ${CASE_DIR}/${ALLRUN_LOGFILE})"
    exit 1
fi

force_time=$(latest_numeric_time "${force_file}" || true)

if [[ -z "${force_time}" ]]; then
    echo "FAIL: the case did not complete in this environment:"
    echo "      the force output contains no time data and the tutorial did not"
    echo "      declare a skip (see ${CASE_DIR}/${ALLRUN_LOGFILE})"
    exit 1
fi

if ! awk "BEGIN {exit !(${force_time} + 0 >= ${REG_END_TIME})}"; then
    echo "FAIL: the force history stops at t = ${force_time}, short of the"
    echo "      requested end time ${REG_END_TIME}: the case did not complete"
    exit 1
fi

tip_uy=$(extract_final_tip_uy)
force_components=$(extract_final_force_components "${force_file}")

if [[ -z "${tip_uy}" ]]; then
    echo "FAIL: Could not extract tip displacement"
    exit 1
fi

if [[ -z "${force_components}" ]]; then
    echo "FAIL: Could not extract the force components from ${force_file}"
    exit 1
fi

final_fx=$(awk '{print $1}' <<< "${force_components}")
final_fy=$(awk '{print $2}' <<< "${force_components}")

tip_uy_diff_abs=$(abs "$(awk "BEGIN {print ${tip_uy} - ${REF_TIP_UY}}")")
final_fx_diff_abs=$(abs "$(awk "BEGIN {print ${final_fx} - ${REF_FX}}")")
final_fy_diff_abs=$(abs "$(awk "BEGIN {print ${final_fy} - ${REF_FY}}")")

failures=0

if awk "BEGIN {exit !(${tip_uy_diff_abs} < ${DISP_TOL})}"; then
    printf "PASS: final tip Uy = %.6g\n" "${tip_uy}"
else
    printf "FAIL: final tip Uy = %.6g\n" "${tip_uy}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_fx_diff_abs} < ${FX_TOL})}"; then
    printf "PASS: final Fx = %.6g\n" "${final_fx}"
else
    printf "FAIL: final Fx = %.6g\n" "${final_fx}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_fy_diff_abs} < ${FY_TOL})}"; then
    printf "PASS: final Fy = %.6g\n" "${final_fy}"
else
    printf "FAIL: final Fy = %.6g\n" "${final_fy}"
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
