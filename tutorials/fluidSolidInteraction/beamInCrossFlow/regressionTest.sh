#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Beam-in-cross-flow FSI regression test
# ============================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

DISP_MAX_TOL=2e-3
DISP_FINAL_TOL=2e-3
FORCE_FINAL_TOL=0.2

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_displacement.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

# Shortened end time for regression efficiency
REGRESSION_END_TIME=1

echo "============================================================"
echo "Beam-in-cross-flow FSI regression test"
echo "Regression end time          = ${REGRESSION_END_TIME}"
echo "Max tip Dx difference       < ${DISP_MAX_TOL}"
echo "Final tip Dx difference     < ${DISP_FINAL_TOL}"
echo "Final total force Fx diff   < ${FORCE_FINAL_TOL}"
echo "============================================================"
echo

extract_max_displacement() {
    awk -v fieldName="Dx" '
    /^#/ {
        for (i = 1; i <= NF; i++) {
            if ($i == fieldName) {
                col = i - 1
            }
        }
        next
    }
    ($1 + 0) == $1 {
        if (col <= 0) {
            col = 2
        }
        if (!seen || $col > maxVal) {
            maxVal = $col
            seen = 1
        }
    }
    END {
        if (seen) print maxVal
    }' "${1}/${DISP_FILE}"
}

extract_final_displacement() {
    awk -v fieldName="Dx" '
    /^#/ {
        for (i = 1; i <= NF; i++) {
            if ($i == fieldName) {
                col = i - 1
            }
        }
        next
    }
    ($1 + 0) == $1 {
        if (col <= 0) {
            col = 3
        }
        last = $col
    }
    END {
        if (last != "") print last
    }' "${1}/${DISP_FILE}"
}

extract_final_force() {
    awk '
    ($1 + 0) == $1 {
        gsub(/[()]/, "", $0)
        if (NF >= 7) {
            last = $2 + $5
        } else {
            last = $2
        }
    }
    END {
        if (last != "") print last
    }' "${1}/${FORCE_FILE}"
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

patch_end_time() {
    local control_dict="$1"
    local tmp_file="${control_dict}.tmp"

    awk -v endTime="${REGRESSION_END_TIME}" '
    /^[[:space:]]*endTime[[:space:]]+/ {
        print "endTime         " endTime ";"
        next
    }
    { print }
    ' "${control_dict}" > "${tmp_file}"

    mv "${tmp_file}" "${control_dict}"
}

prepare_case() {
    local coupling="$1"
    local case_dir="${REGRESSION_ROOT}/${coupling}"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        local base_item
        base_item=$(basename "${item}")

        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi

        cp -a "${item}" "${case_dir}/"
    done

    (
        cd "${case_dir}"

        ln -vnsf "fsiProperties.${coupling}" constant/fsiProperties
        ln -vnsf "solidProperties.${coupling}" constant/solid/solidProperties
        ln -vnsf "controlDict.${coupling}" system/controlDict
        ln -vnsf "fvSolution.${coupling}" system/fluid/fvSolution

        patch_end_time "system/controlDict.${coupling}"
    ) > /dev/null

    echo "${case_dir}"
}

run_case() {
    local case_dir="$1"
    local coupling="$2"

    echo "Running ${coupling} regression case"

    (
        cd "${case_dir}"
        ./Allclean > /dev/null 2>&1 || true
        ./Allrun "${coupling}" > "${ALLRUN_LOGFILE}" 2>&1

        mkdir -p postProcessing/fluid/forces/0
        cd postProcessing/fluid/forces/0

        if [[ ! -e force.dat && -f ../../../../forces/0/forces.dat ]]; then
            ln -s ../../../../forces/0/forces.dat force.dat
        fi

        if [[ ! -e force.dat && -f forces.dat ]]; then
            ln -s forces.dat force.dat
        fi
    )
}

check_case() {
    local coupling="$1"
    local case_dir="$2"
    local ref_max_disp="$3"
    local ref_final_disp="$4"
    local ref_final_force="$5"

    local max_disp
    local final_disp
    local final_force
    local disp_diff
    local disp_diff_abs
    local final_disp_diff
    local final_disp_diff_abs
    local final_force_diff
    local final_force_diff_abs
    local failures=0

    max_disp=$(extract_max_displacement "${case_dir}")
    final_disp=$(extract_final_displacement "${case_dir}")
    final_force=$(extract_final_force "${case_dir}")

    if [[ -z "${max_disp}" || -z "${final_disp}" || -z "${final_force}" ]]; then
        echo "FAIL [${coupling}]: Could not extract regression quantities"
        return 1
    fi

    disp_diff=$(awk "BEGIN {print ${max_disp} - ${ref_max_disp}}")
    disp_diff_abs=$(abs "${disp_diff}")

    final_disp_diff=$(awk "BEGIN {print ${final_disp} - ${ref_final_disp}}")
    final_disp_diff_abs=$(abs "${final_disp_diff}")

    final_force_diff=$(awk "BEGIN {print ${final_force} - ${ref_final_force}}")
    final_force_diff_abs=$(abs "${final_force_diff}")

    if awk "BEGIN {exit !(${disp_diff_abs} < ${DISP_MAX_TOL})}"; then
        printf "PASS [%s]: max displacement = %.6g (Δ = %.3g)\n" \
            "${coupling}" "${max_disp}" "${disp_diff_abs}"
    else
        printf "FAIL [%s]: max displacement = %.6g (Δ = %.3g)\n" \
            "${coupling}" "${max_disp}" "${disp_diff_abs}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${final_disp_diff_abs} < ${DISP_FINAL_TOL})}"; then
        printf "PASS [%s]: final tip Dx = %.6g (Δ = %.3g)\n" \
            "${coupling}" "${final_disp}" "${final_disp_diff_abs}"
    else
        printf "FAIL [%s]: final tip Dx = %.6g (Δ = %.3g)\n" \
            "${coupling}" "${final_disp}" "${final_disp_diff_abs}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${final_force_diff_abs} < ${FORCE_FINAL_TOL})}"; then
        printf "PASS [%s]: final total Fx = %.6g (Δ = %.3g)\n" \
            "${coupling}" "${final_force}" "${final_force_diff_abs}"
    else
        printf "FAIL [%s]: final total Fx = %.6g (Δ = %.3g)\n" \
            "${coupling}" "${final_force}" "${final_force_diff_abs}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

REF_MAX_DISP=0.039909
REF_FINAL_DISP=0.039232
REF_FINAL_FORCE=6.54692

aitken_case=$(prepare_case aitken)
run_case "${aitken_case}" aitken

iqnils_case=$(prepare_case iqnils)
run_case "${iqnils_case}" iqnils

echo

failures=0

check_case aitken "${aitken_case}" \
    "${REF_MAX_DISP}" "${REF_FINAL_DISP}" "${REF_FINAL_FORCE}" \
    || failures=$((failures + $?))

check_case iqnils "${iqnils_case}" \
    "${REF_MAX_DISP}" "${REF_FINAL_DISP}" "${REF_FINAL_FORCE}" \
    || failures=$((failures + $?))

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
