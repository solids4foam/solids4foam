#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# ============================================================
# ringAddedMass FSI regression test
# ============================================================

# Shortened regression horizon: about 1.05 wet periods (104 time-steps), which
# contains two zero crossings of the ovalling displacement
REG_END_TIME=3.64

# Exact half period of the wet ring, pi/omega, for the tutorial's fluid density
# (see verification/reference/ringAddedMass_verification_references.json)
EXACT_HALF_PERIOD=1.740964

# Regression tolerances
HALF_PERIOD_EXACT_TOL=0.01    # relative difference from the exact half period
HALF_PERIOD_TOL=2e-4          # relative difference from the reference
DISP_TOL=2e-6                 # final displacement absolute difference (m)

# Reference values at REG_END_TIME (OpenFOAM-v2512, with foam-extend values
# selected below)
ref_half_period() {
    if [[ "${WM_PROJECT:-}" == "foam" ]]; then
        case "$1" in
            iqnils) echo 1.7411818 ;;
            robin) echo 1.7436173 ;;
        esac
        return
    fi
    case "$1" in
        iqnils) echo 1.7420511 ;;
        robin) echo 1.7432661 ;;
    esac
}

ref_final_ux() {
    if [[ "${WM_PROJECT:-}" == "foam" ]]; then
        case "$1" in
            iqnils) echo 0.000155080 ;;
            robin) echo 0.000112373 ;;
        esac
        return
    fi
    case "$1" in
        iqnils) echo 0.000121865 ;;
        robin) echo 0.000120717 ;;
    esac
}

VARIANTS=(iqnils robin)

ALLRUN_LOGFILE="log.Allrun"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDispX.dat"

echo "============================================================"
echo "ringAddedMass FSI regression test"
echo "Regression end time               = ${REG_END_TIME}"
echo "Half period vs exact tolerance    < ${HALF_PERIOD_EXACT_TOL}"
echo "Half period vs reference tolerance < ${HALF_PERIOD_TOL}"
echo "Final displacement tolerance      < ${DISP_TOL}"
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

final_ux() {
    awk '($1 + 0) == $1 { ux = $2 } END { if (ux != "") print ux }' "$1"
}

# Difference between the first two zero-crossing times of u_x, i.e. half a
# period, from linear interpolation between the samples. The ring starts from
# its undeformed shape, so the first crossing is after half a period
half_period() {
    awk '
    ($1 + 0) == $1 {
        if (n > 0 && ux * $2 < 0) {
            crossing[++k] = t + ux / (ux - $2) * ($1 - t)
        }
        t = $1
        ux = $2
        n++
    }
    END {
        if (k >= 2) printf "%.8g\n", crossing[2] - crossing[1]
    }' "$1"
}

# Ensure GNU sed is available and resolved into SOLIDS4FOAM_SED
solids4Foam::requireGnuSed

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

    half=$(half_period "${disp_file}")
    ux=$(final_ux "${disp_file}")
    if [[ -z "${half}" || -z "${ux}" ]]; then
        echo "FAIL (${variant}): could not extract the half period"
        failures=$((failures + 1))
        continue
    fi

    ref_half=$(ref_half_period "${variant}")
    ref_ux=$(ref_final_ux "${variant}")
    exact_diff=$(abs "$(awk "BEGIN {print ${half}/${EXACT_HALF_PERIOD} - 1}")")
    ref_diff=$(abs "$(awk "BEGIN {print ${half}/${ref_half} - 1}")")
    ux_diff=$(abs "$(awk "BEGIN {print ${ux} - ${ref_ux}}")")

    if awk "BEGIN {exit !(${exact_diff} < ${HALF_PERIOD_EXACT_TOL})}"; then
        printf "PASS (%s): half period = %.6g s, exact %.6g s (%.2g)\n" \
            "${variant}" "${half}" "${EXACT_HALF_PERIOD}" "${exact_diff}"
    else
        printf "FAIL (%s): half period = %.6g s, exact %.6g s (%.2g)\n" \
            "${variant}" "${half}" "${EXACT_HALF_PERIOD}" "${exact_diff}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${ref_diff} < ${HALF_PERIOD_TOL})}"; then
        printf "PASS (%s): half period matches the reference (%.2g)\n" \
            "${variant}" "${ref_diff}"
    else
        printf "FAIL (%s): half period = %.8g s, reference %.8g s (%.2g)\n" \
            "${variant}" "${half}" "${ref_half}" "${ref_diff}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${ux_diff} < ${DISP_TOL})}"; then
        printf "PASS (%s): final u_x = %.6g m (diff %.2g)\n" \
            "${variant}" "${ux}" "${ux_diff}"
    else
        printf "FAIL (%s): final u_x = %.6g m, reference %.6g m (diff %.2g)\n" \
            "${variant}" "${ux}" "${ref_ux}" "${ux_diff}"
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
