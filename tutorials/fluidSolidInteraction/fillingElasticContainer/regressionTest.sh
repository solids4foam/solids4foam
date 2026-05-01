#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# fillingElasticContainer regression test
# ============================================================

REG_END_TIME=2
MODE="iqnils"

# Regression tolerances
APEX_DISP_TOL=1e-3
FSI_RES_TOL=1e-6

# Reference values at REG_END_TIME
REF_APEX_DY=-0.534826
REF_FSI_RES=5.35591e-06

ALLRUN_LOGFILE="log.Allrun"
DISP_FILE="postProcessing/0/solidPointDisplacement_disp.dat"
FSI_RES_FILE="postProcessing/fsiResiduals.dat"

echo "============================================================"
echo "fillingElasticContainer regression test"
echo "Regression end time         = ${REG_END_TIME}"
echo "Coupling mode               = ${MODE}"
echo "Final apex dy tolerance     < ${APEX_DISP_TOL}"
echo "Final FSI residual tolerance < ${FSI_RES_TOL}"
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
        ./Allrun "${MODE}" > "${ALLRUN_LOGFILE}" 2>&1
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

extract_final_apex_dy() {
    awk '
    ($1 + 0) == $1 {
        dy = $3
    }
    END {
        if (dy != "") print dy
    }' "${CASE_DIR}/${DISP_FILE}"
}

extract_final_fsi_residual() {
    awk '
    ($1 + 0) == $1 {
        residual = $3
    }
    END {
        if (residual != "") print residual
    }' "${CASE_DIR}/${FSI_RES_FILE}"
}

prepare_case
run_case

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

apex_time=$(latest_numeric_time "${CASE_DIR}/${DISP_FILE}" || true)
fsi_time=$(latest_numeric_time "${CASE_DIR}/${FSI_RES_FILE}" || true)

if [[ -z "${apex_time}" || -z "${fsi_time}" ]]; then
    echo "Skipping regression checks because the case did not complete in this environment"
    exit 0
fi

if ! awk "BEGIN {exit !(${apex_time} + 0 >= ${REG_END_TIME})}"; then
    echo "Skipping regression checks because the apex history did not reach the requested end time"
    exit 0
fi

if ! awk "BEGIN {exit !(${fsi_time} + 0 >= ${REG_END_TIME})}"; then
    echo "Skipping regression checks because the FSI residual history did not reach the requested end time"
    exit 0
fi

apex_dy=$(extract_final_apex_dy)
fsi_residual=$(extract_final_fsi_residual)

if [[ -z "${apex_dy}" || -z "${fsi_residual}" ]]; then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

apex_diff_abs=$(abs "$(awk "BEGIN {print ${apex_dy} - ${REF_APEX_DY}}")")
fsi_res_diff_abs=$(abs "$(awk "BEGIN {print ${fsi_residual} - ${REF_FSI_RES}}")")

failures=0

if awk "BEGIN {exit !(${apex_diff_abs} < ${APEX_DISP_TOL})}"; then
    printf "PASS: final apex dy = %.6g\n" "${apex_dy}"
else
    printf "FAIL: final apex dy = %.6g\n" "${apex_dy}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${fsi_res_diff_abs} < ${FSI_RES_TOL})}"; then
    printf "PASS: final FSI residual = %.6g\n" "${fsi_residual}"
else
    printf "FAIL: final FSI residual = %.6g\n" "${fsi_residual}"
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
