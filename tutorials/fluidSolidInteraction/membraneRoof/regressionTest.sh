#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# membraneRoof FSI regression test
# ============================================================

# Shortened regression horizon: the first 10 time-steps of the inlet ramp,
# during which the roof response is smooth. A full run takes about 18 min on
# 6 processes.
REG_END_TIME=0.2

# Regression tolerances
DISP_TOL=1e-4
FY_TOL=20

# Reference values at REG_END_TIME: vertical displacement of the roof centre
# and the total vertical force on the fluid side of the roof, from a serial
# OpenFOAM-v2412 run on Linux. A 6-process OpenFOAM-v2412 run on macOS differs
# by 3.8e-6 m and 1.4 N.
REF_UY=-0.0610105
REF_FY=-28041.0

ALLRUN_LOGFILE="log.Allrun"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

echo "============================================================"
echo "membraneRoof FSI regression test"
echo "Regression end time         = ${REG_END_TIME}"
echo "Roof centre Uy tolerance    < ${DISP_TOL}"
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

    # Set the end time with awk rather than sed -i, which differs between GNU
    # and BSD sed
    local controlDict="${CASE_DIR}/system/controlDict"
    awk -v endTime="${REG_END_TIME}" '
        /^endTime[[:space:]]/ { $0 = "endTime         " endTime ";" }
        { print }
    ' "${controlDict}" > "${controlDict}.tmp"
    mv "${controlDict}.tmp" "${controlDict}"
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

extract_final_uy() {
    awk '
    ($1 + 0) == $1 {
        uy = $3
    }
    END {
        if (uy != "") print uy
    }' "${CASE_DIR}/${DISP_FILE}"
}

extract_final_fy() {
    # OpenFOAM.com writes the time and the total force first
    awk '
    ($1 + 0) == $1 {
        gsub(/[()]/, "", $0)
        fy = $3
    }
    END {
        if (fy != "") print fy
    }' "${CASE_DIR}/${FORCE_FILE}"
}

prepare_case
run_case

# A skip is only valid if the tutorial declared one in the Allrun log. Anything
# else that leaves the expected output missing or incomplete is a failure.
if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

for file in "${DISP_FILE}" "${FORCE_FILE}"; do
    if [[ ! -f "${CASE_DIR}/${file}" ]]; then
        echo "FAIL: ${file} is missing and the tutorial did not declare a skip"
        echo "      (see ${CASE_DIR}/${ALLRUN_LOGFILE})"
        exit 1
    fi

    last_time=$(latest_numeric_time "${CASE_DIR}/${file}" || true)

    if [[ -z "${last_time}" ]] || \
        ! awk "BEGIN {exit !(${last_time} + 0 >= ${REG_END_TIME})}"; then
        echo "FAIL: ${file} stops at t = ${last_time:-none}, short of the"
        echo "      requested end time ${REG_END_TIME}: the case did not complete"
        exit 1
    fi
done

final_uy=$(extract_final_uy)
final_fy=$(extract_final_fy)

if [[ -z "${final_uy}" || -z "${final_fy}" ]]; then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

uy_diff_abs=$(abs "$(awk "BEGIN {print ${final_uy} - ${REF_UY}}")")
fy_diff_abs=$(abs "$(awk "BEGIN {print ${final_fy} - ${REF_FY}}")")

failures=0

if awk "BEGIN {exit !(${uy_diff_abs} < ${DISP_TOL})}"; then
    printf "PASS: final roof centre Uy = %.6g (Δ = %.3g)\n" "${final_uy}" "${uy_diff_abs}"
else
    printf "FAIL: final roof centre Uy = %.6g (Δ = %.3g)\n" "${final_uy}" "${uy_diff_abs}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${fy_diff_abs} < ${FY_TOL})}"; then
    printf "PASS: final Fy = %.6g (Δ = %.3g)\n" "${final_fy}" "${fy_diff_abs}"
else
    printf "FAIL: final Fy = %.6g (Δ = %.3g)\n" "${final_fy}" "${fy_diff_abs}"
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
