#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# membraneRoof FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

REGRESSION_END_TIME=2

REF_FINAL_DY=-0.39024
REF_FINAL_FORCE_Y=-204924.3

DY_TOL=0.01
FORCE_Y_TOL=5000

ALLRUN_LOGFILE="log.Allrun"
ALLWMAKE_LOGFILE="log.Allwmake"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"

echo "============================================================"
echo "membraneRoof FSI regression test"
echo "Regression end time          = ${REGRESSION_END_TIME}"
echo "Final point Dy difference    < ${DY_TOL}"
echo "Final membrane Fy difference < ${FORCE_Y_TOL}"
echo "============================================================"
echo

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
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        local base_item
        base_item=$(basename "${item}")

        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi

        cp -a "${item}" "${CASE_DIR}/"
    done

    patch_end_time "${CASE_DIR}/system/controlDict"
}

find_force_file() {
    local case_dir="$1"
    local candidate

    for candidate in \
        "${case_dir}/postProcessing/fluid/forces/0/force.dat" \
        "${case_dir}/postProcessing/fluid/forces/0/forces.dat" \
        "${case_dir}/postProcessing/forces/0/force.dat" \
        "${case_dir}/postProcessing/forces/0/forces.dat"
    do
        if [[ -f "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done

    return 1
}

extract_final_dy() {
    local case_dir="$1"

    awk -v fieldName="Dy" '
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
        val = $col
    }
    END {
        if (val != "") print val
    }' "${case_dir}/${DISP_FILE}"
}

extract_final_force_y() {
    local case_dir="$1"
    local force_file

    force_file=$(find_force_file "${case_dir}") || return 1

    awk -v fieldName="total_y" '
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
        val = $col
    }
    END {
        if (val != "") print val
    }' "${force_file}"
}

check_results() {
    local case_dir="$1"
    local failures=0

    if [[ ! -f "${case_dir}/${DISP_FILE}" ]]; then
        echo "FAIL: Could not find ${DISP_FILE}"
        return 1
    fi

    local final_dy
    final_dy=$(extract_final_dy "${case_dir}")

    local final_force_y
    final_force_y=$(extract_final_force_y "${case_dir}")

    if [[ -z "${final_dy}" || -z "${final_force_y}" ]]; then
        echo "FAIL: Could not extract regression quantities"
        return 1
    fi

    local dy_diff
    dy_diff=$(awk "BEGIN {print ${final_dy} - ${REF_FINAL_DY}}")
    local dy_diff_abs
    dy_diff_abs=$(abs "${dy_diff}")

    local force_y_diff
    force_y_diff=$(awk "BEGIN {print ${final_force_y} - ${REF_FINAL_FORCE_Y}}")
    local force_y_diff_abs
    force_y_diff_abs=$(abs "${force_y_diff}")

    if awk "BEGIN {exit !(${dy_diff_abs} < ${DY_TOL})}"; then
        printf "PASS: final Dy = %.6g (Δ = %.3g)\n" "${final_dy}" "${dy_diff_abs}"
    else
        printf "FAIL: final Dy = %.6g (Δ = %.3g)\n" "${final_dy}" "${dy_diff_abs}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${force_y_diff_abs} < ${FORCE_Y_TOL})}"; then
        printf "PASS: final Fy = %.6g (Δ = %.3g)\n" \
            "${final_force_y}" "${force_y_diff_abs}"
    else
        printf "FAIL: final Fy = %.6g (Δ = %.3g)\n" \
            "${final_force_y}" "${force_y_diff_abs}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

if [[ "${CHECK_ONLY}" == false ]]; then
    prepare_case

    run_status=0

    set +e
    (
        cd "${CASE_DIR}"
        ./Allclean > /dev/null 2>&1 || true
        (
            cd src
            ./Allwmake > "../${ALLWMAKE_LOGFILE}" 2>&1
        )
        ./Allrun > "${ALLRUN_LOGFILE}" 2>&1
    )
    run_status=$?
    set -e

    if declare -F "solids4Foam::regressionCaseSkipped" > /dev/null \
        && solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"
    then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    if (( run_status != 0 )); then
        echo "FAIL: ./Allrun failed in ${CASE_DIR}"
        echo "Check ${CASE_DIR}/${ALLWMAKE_LOGFILE}, ${CASE_DIR}/${ALLRUN_LOGFILE}, and ${CASE_DIR}/log.solids4Foam"
        exit 1
    fi
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
    CASE_DIR="${SCRIPT_DIR}"
fi

failures=0
check_results "${CASE_DIR}" || failures=$((failures + $?))

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
