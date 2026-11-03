#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# ============================================================
# radialTorsionalRing FSI regression test
# ============================================================

# Shortened regression horizon: 120 time-steps of the M1 mesh, combined
# drive (radial and torsional), serial. Only the start of the drive is run, so
# this is a check that the tutorial launches and reproduces the same discrete
# solution, not a verification of the benchmark (see verification/README.md).
REG_END_TIME=0.3

# Regression tolerances. The interface radial displacement u_B at t = 0.3 s is
# 1.0928e-3 m in the solver and 1.1008e-3 m in the semi-analytical reference
# (an M1-mesh difference of 0.7%, 8.0e-6 m). The value reproduces the
# accepted M1 combined run of the study to round-off.
EXACT_TOL=3e-5          # |u_B - u_B(semi-analytical reference)| (m)
DISP_TOL=1e-7           # |u_B - u_B(regression reference)| (m)
ALPHA_TOL=2e-8          # |alpha_B - alpha_B(regression reference)| (rad)

# Reference values at REG_END_TIME (OpenFOAM-v2512, M1, combined,
# nNonOrthogonalCorrectors 2, 1 rank)
REF_U_B=0.001092801156254323
REF_ALPHA_B=-1.218218794956131e-05

ALLRUN_LOGFILE="log.Allrun"
DIAG_FILE="rtDiagSolid.dat"

echo "============================================================"
echo "radialTorsionalRing FSI regression test"
echo "Regression end time                 = ${REG_END_TIME}"
echo "u_B vs semi-analytical reference    < ${EXACT_TOL}"
echo "u_B vs regression reference         < ${DISP_TOL}"
echo "alpha_B vs regression reference     < ${ALPHA_TOL}"
echo "============================================================"
echo

prepare_case() {
    local case_dir="$1"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        if [[ "${base_item}" == "verification" ]]; then
            # only the generator and analysis scripts, not any study cases
            mkdir -p "${case_dir}/verification"
            cp -a "${item}/scripts" "${case_dir}/verification/"
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

failures=0
case_dir="${REGRESSION_ROOT}/m1_combined"
echo "Running M1 combined"
prepare_case "${case_dir}"
(
    cd "${case_dir}"
    ./Allclean > /dev/null 2>&1 || true
    ./Allrun M1 combined serial "endTime=${REG_END_TIME}" > "${ALLRUN_LOGFILE}" 2>&1
)

# A skip is only valid if the tutorial declared one in the Allrun log
if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

solver_log="${case_dir}/log.solids4Foam"
if [[ ! -f "${solver_log}" ]] \
    || grep -q "FOAM FATAL" "${solver_log}" \
    || ! grep -q "^End" "${solver_log}"; then
    echo "FAIL: the solver did not run to completion"
    echo "      (see ${case_dir}/${ALLRUN_LOGFILE})"
    exit 1
fi

diag_file="${case_dir}/${DIAG_FILE}"
if [[ ! -f "${diag_file}" ]]; then
    echo "FAIL: the diagnostics output ${DIAG_FILE} is missing"
    exit 1
fi

# Last row of rtDiagSolid.dat: time, u_B (column 3), alpha_B (column 4); the
# semi-analytical reference at the same time and face radius
IFS=" " read -r end_time ur alpha < <(awk 'NF >= 4 { t=$1; u=$3; a=$4 } END { print t, u, a }' "${diag_file}")
if ! awk "BEGIN {exit !(${end_time} + 1e-9 >= ${REG_END_TIME})}"; then
    echo "FAIL: the diagnostics stop at t = ${end_time}, short of ${REG_END_TIME}"
    exit 1
fi
ref_ur=$(python3 -B -I -c "
import json, sys
r = json.load(open('${case_dir}/reference.json'))
i = min(range(len(r['t'])), key=lambda k: abs(r['t'][k] - ${end_time}))
print('%.12g' % r['u_Rf'][i])")

exact_diff=$(abs "$(awk "BEGIN {print ${ur} - ${ref_ur}}")")
ref_diff=$(abs "$(awk "BEGIN {print ${ur} - ${REF_U_B}}")")
alpha_diff=$(abs "$(awk "BEGIN {print ${alpha} - (${REF_ALPHA_B})}")")

check() {
    # name, difference, tolerance
    if awk "BEGIN {exit !($2 < $3)}"; then
        printf "PASS: %s (diff %.2g < %s)\n" "$1" "$2" "$3"
    else
        printf "FAIL: %s (diff %.2g, tolerance %s)\n" "$1" "$2" "$3"
        failures=$((failures + 1))
    fi
}

printf "t = %s, u_B = %.10g m, alpha_B = %.6g rad\n" "${end_time}" "${ur}" "${alpha}"
check "u_B against the semi-analytical reference" "${exact_diff}" "${EXACT_TOL}"
check "u_B against the regression reference" "${ref_diff}" "${DISP_TOL}"
check "alpha_B against the regression reference" "${alpha_diff}" "${ALPHA_TOL}"
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
