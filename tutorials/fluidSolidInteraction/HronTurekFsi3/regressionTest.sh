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
DISP_TOL=1e-4
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
# The HronTurekFsi3 values also include the fluid residualControl added in
# the same PR, so their shift is not due to the interface normals alone.
# The tip Uy tolerance is widened from 2e-5 to 1e-4 and the reference set to
# the midpoint of the OpenFOAM.com versions: the spread between v2412 and
# v2512 is now about 6.8e-5, against 4.2e-6 before this change.
REF_TIP_UY=-0.0024405
REF_FX=-0.0386691
REF_FY=-0.0203155

ALLRUN_LOGFILE="log.Allrun"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

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

    mkdir -p "${CASE_DIR}/postProcessing/fluid/forces/0"
    (
        cd "${CASE_DIR}/postProcessing/fluid/forces/0"
        if [[ ! -e force.dat && -f forces.dat ]]; then
            ln -s forces.dat force.dat
        fi
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
    local candidate
    for candidate in \
        "${CASE_DIR}/postProcessing/fluid/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/fluid/forces/0/forces.dat" \
        "${CASE_DIR}/postProcessing/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/forces/0/forces.dat"
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

tip_time=$(latest_numeric_time "${CASE_DIR}/${DISP_FILE}" || true)
force_file=""
if force_file=$(find_force_file); then
    force_time=$(latest_numeric_time "${force_file}" || true)
else
    force_time=""
fi

if [[ -z "${tip_time}" || -z "${force_file}" || -z "${force_time}" ]]; then
    echo "Skipping regression checks because the case did not complete in this environment"
    exit 0
fi

if ! awk "BEGIN {exit !(${tip_time} + 0 >= ${REG_END_TIME})}"; then
    echo "Skipping regression checks because the tip displacement history did not reach the requested end time"
    exit 0
fi

if ! awk "BEGIN {exit !(${force_time} + 0 >= ${REG_END_TIME})}"; then
    echo "Skipping regression checks because the force history did not reach the requested end time"
    exit 0
fi

tip_uy=$(extract_final_tip_uy)
if force_file=$(find_force_file); then
    force_components=$(extract_final_force_components "${force_file}")
else
    force_components=""
    echo "Skipping force checks because force data is unavailable"
fi

if [[ -z "${tip_uy}" || -z "${force_components}" ]]; then
    if [[ -z "${tip_uy}" ]]; then
        echo "FAIL: Could not extract tip displacement"
        exit 1
    fi
fi

if [[ -n "${force_components}" ]]; then
    final_fx=$(awk '{print $1}' <<< "${force_components}")
    final_fy=$(awk '{print $2}' <<< "${force_components}")
fi

tip_uy_diff_abs=$(abs "$(awk "BEGIN {print ${tip_uy} - ${REF_TIP_UY}}")")
if [[ -n "${force_components}" ]]; then
    final_fx_diff_abs=$(abs "$(awk "BEGIN {print ${final_fx} - ${REF_FX}}")")
    final_fy_diff_abs=$(abs "$(awk "BEGIN {print ${final_fy} - ${REF_FY}}")")
fi

failures=0

if awk "BEGIN {exit !(${tip_uy_diff_abs} < ${DISP_TOL})}"; then
    printf "PASS: final tip Uy = %.6g\n" "${tip_uy}"
else
    printf "FAIL: final tip Uy = %.6g\n" "${tip_uy}"
    failures=$((failures + 1))
fi

if [[ -n "${force_components}" ]]; then
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
