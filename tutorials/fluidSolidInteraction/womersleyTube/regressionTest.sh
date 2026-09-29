#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# ============================================================
# womersleyTube FSI regression test
# ============================================================

# Shortened regression horizon: half a period (50 time-steps)
REG_END_TIME=25

# Exact radial wall displacement at x = L/2, Re[a exp(i omega t)], with
# |a| and arg(a) from verification/reference/womersleyTube_verification_references.json
EXACT_AMPLITUDE=4.512494768325512e-4
EXACT_PHASE=-1.0453187730902433
OMEGA=0.12566370614359174

# Regression tolerances. The difference from the exact solution, 2.6e-6 m
# (IQN-ILS) and 4.6e-6 m (Robin) of an amplitude of 4.5e-4 m, is mostly the
# start-up transient, which decays over the first periods
EXACT_TOL=1e-5          # final displacement difference from the exact (m)
DISP_TOL=1e-6           # final displacement difference from the reference (m)

# Reference values at REG_END_TIME (OpenFOAM-v2412)
ref_final_ur() {
    case "$1" in
        iqnils) echo -0.0002237665584 ;;
        robin) echo -0.0002217864674 ;;
    esac
}

VARIANTS=(iqnils robin)

ALLRUN_LOGFILE="log.Allrun"
DISP_FILE="postProcessing/0/solidPointDisplacement_wallMid.dat"

echo "============================================================"
echo "womersleyTube FSI regression test"
echo "Regression end time              = ${REG_END_TIME}"
echo "Displacement vs exact tolerance  < ${EXACT_TOL}"
echo "Displacement vs reference tolerance < ${DISP_TOL}"
echo "============================================================"
echo

prepare_case() {
    local case_dir="$1"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" \
           || "${base_item}" == "verification" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done

    "${SOLIDS4FOAM_SED}" -i \
        "s/^\(endTime[[:space:]]*\).*/\1${REG_END_TIME};/" \
        "${case_dir}/system/controlDict"
}

run_case() {
    local case_dir="$1"
    local variant="$2"
    (
        cd "${case_dir}"
        ./Allclean > /dev/null 2>&1 || true
        ./Allrun "${variant}" > "${ALLRUN_LOGFILE}" 2>&1
    )
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

latest_numeric_time() {
    awk '($1 + 0) == $1 { time = $1 } END { if (time != "") print time }' "$1"
}

# Radial displacement of the last sample; the point is on the wedge plane at
# -0.5 degrees
final_ur() {
    awk '
    ($1 + 0) == $1 {
        th = -0.5*atan2(0, -1)/180
        ur = $3*cos(th) + $4*sin(th)
    }
    END { if (ur != "") printf "%.10g\n", ur }' "$1"
}

# Ensure GNU sed is available and resolved into SOLIDS4FOAM_SED
solids4Foam::requireGnuSed

exact_ur=$(awk -v a="${EXACT_AMPLITUDE}" -v p="${EXACT_PHASE}" \
    -v w="${OMEGA}" -v t="${REG_END_TIME}" 'BEGIN {printf "%.10g", a*cos(w*t + p)}')

failures=0

for variant in "${VARIANTS[@]}"; do
    case_dir="${REGRESSION_ROOT}/${variant}"
    echo "Running ${variant}"
    prepare_case "${case_dir}"
    run_case "${case_dir}" "${variant}"

    # A skip is only valid if the tutorial declared one in the Allrun log
    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    solver_log="${case_dir}/log.solids4Foam"
    if [[ ! -f "${solver_log}" ]] \
        || grep -q "FOAM FATAL" "${solver_log}" \
        || ! grep -q "^End" "${solver_log}"; then
        echo "FAIL (${variant}): the solver did not run to completion"
        echo "      (see ${case_dir}/${ALLRUN_LOGFILE})"
        failures=$((failures + 1))
        continue
    fi

    disp_file="${case_dir}/${DISP_FILE}"
    if [[ ! -f "${disp_file}" ]]; then
        echo "FAIL (${variant}): displacement output is missing"
        echo "      (see ${case_dir}/${ALLRUN_LOGFILE})"
        failures=$((failures + 1))
        continue
    fi

    end_time=$(latest_numeric_time "${disp_file}" || true)
    if [[ -z "${end_time}" ]] \
        || ! awk "BEGIN {exit !(${end_time} + 1e-9 >= ${REG_END_TIME})}"; then
        echo "FAIL (${variant}): the displacement history stops at"
        echo "      t = ${end_time:-none}, short of ${REG_END_TIME}"
        failures=$((failures + 1))
        continue
    fi

    ur=$(final_ur "${disp_file}")
    if [[ -z "${ur}" ]]; then
        echo "FAIL (${variant}): could not extract the wall displacement"
        failures=$((failures + 1))
        continue
    fi

    ref_ur=$(ref_final_ur "${variant}")
    exact_diff=$(abs "$(awk "BEGIN {print ${ur} - ${exact_ur}}")")
    ref_diff=$(abs "$(awk "BEGIN {print ${ur} - ${ref_ur}}")")

    if awk "BEGIN {exit !(${exact_diff} < ${EXACT_TOL})}"; then
        printf "PASS (%s): final u_r = %.6g m, exact %.6g m (diff %.2g)\n" \
            "${variant}" "${ur}" "${exact_ur}" "${exact_diff}"
    else
        printf "FAIL (%s): final u_r = %.6g m, exact %.6g m (diff %.2g)\n" \
            "${variant}" "${ur}" "${exact_ur}" "${exact_diff}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${ref_diff} < ${DISP_TOL})}"; then
        printf "PASS (%s): final u_r matches the reference (diff %.2g)\n" \
            "${variant}" "${ref_diff}"
    else
        printf "FAIL (%s): final u_r = %.8g m, reference %.8g m (diff %.2g)\n" \
            "${variant}" "${ur}" "${ref_ur}" "${ref_diff}"
        failures=$((failures + 1))
    fi
    echo
done

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
