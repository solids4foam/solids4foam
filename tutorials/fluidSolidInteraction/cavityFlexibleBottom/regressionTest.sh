#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# ============================================================
# cavityFlexibleBottom FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression settings (shared by partitioned and monolithic)
# ------------------------------------------------------------

# Shortened end time: both methods are near steady state by t=100
# (partitioned Dy=-0.218, monolithic Dy=-0.212 at t=100)
REGRESSION_END_TIME=100

# Number of samples from end of force.dat to average (covers t=90.5-100)
FORCE_AVG_SAMPLES=20

# Reference values: near-steady-state values common to both methods
REF_FINAL_DY=-0.22
REF_MEAN_FORCE=-5150

# Tolerances: loosened slightly to accommodate both methods at t=100
DISP_TOL=0.010       # covers partitioned (-0.218) and monolithic (-0.212)
FORCE_TOL=150        # covers partitioned (-5150) and monolithic (-5117)

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

# ------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------

CHECK_ONLY=false
RUN_MONOLITHIC=false

for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run)
            CHECK_ONLY=true
            ;;
        --monolithic)
            RUN_MONOLITHIC=true
            ;;
        *)
            ;;
    esac
done

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

abs() { awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'; }

patch_end_time() {
    local control_dict="$1"
    local new_end_time="$2"
    local tmp_file="${control_dict}.tmp"
    awk -v endTime="${new_end_time}" '
    /^[[:space:]]*endTime[[:space:]]+/ { print "endTime         " endTime ";"; next }
    { print }
    ' "${control_dict}" > "${tmp_file}"
    mv "${tmp_file}" "${control_dict}"
}

extract_vertical_displacement() {
    local case_dir="$1"
    awk '{print $3}' "${case_dir}/${DISP_FILE}" | tail -1
}

extract_mean_force() {
    local case_dir="$1"
    tail -n "${FORCE_AVG_SAMPLES}" "${case_dir}/${FORCE_FILE}" | \
    awk '
    ($1+0)==$1 {
        gsub(/[()]/, "", $0)
        sum += $3
        n++
    }
    END {
        if (n > 0) print sum/n
    }'
}

setup_force_dat() {
    local case_dir="$1"
    mkdir -p "${case_dir}/postProcessing/fluid/forces/0"
    (
        cd "${case_dir}/postProcessing/fluid/forces/0"
        if [[ ! -e force.dat && -f ../../../../forces/0/forces.dat ]]; then
            ln -s ../../../../forces/0/forces.dat force.dat
        fi
        if [[ ! -e force.dat && -f forces.dat ]]; then
            ln -s forces.dat force.dat
        fi
    )
}

check_results() {
    local label="$1"
    local case_dir="$2"
    local failures=0

    local final_dy
    final_dy=$(extract_vertical_displacement "${case_dir}")

    local mean_force
    mean_force=$(extract_mean_force "${case_dir}")

    if [[ -z "${final_dy}" || -z "${mean_force}" ]]; then
        echo "FAIL [${label}]: Could not extract regression quantities"
        return 1
    fi

    local dy_diff
    dy_diff=$(awk "BEGIN {print ${final_dy} - ${REF_FINAL_DY}}")
    local dy_diff_abs
    dy_diff_abs=$(abs "${dy_diff}")

    local force_diff
    force_diff=$(awk "BEGIN {print ${mean_force} - ${REF_MEAN_FORCE}}")
    local force_diff_abs
    force_diff_abs=$(abs "${force_diff}")

    if awk "BEGIN {exit !(${dy_diff_abs} < ${DISP_TOL})}"; then
        printf "PASS [%s]: final Dy = %.6g (Δ = %.3g)\n" \
            "${label}" "${final_dy}" "${dy_diff_abs}"
    else
        printf "FAIL [%s]: final Dy = %.6g (Δ = %.3g)\n" \
            "${label}" "${final_dy}" "${dy_diff_abs}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${force_diff_abs} < ${FORCE_TOL})}"; then
        printf "PASS [%s]: mean Fy = %.6g (Δ = %.3g)\n" \
            "${label}" "${mean_force}" "${force_diff_abs}"
    else
        printf "FAIL [%s]: mean Fy = %.6g (Δ = %.3g)\n" \
            "${label}" "${mean_force}" "${force_diff_abs}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

# ------------------------------------------------------------
# Monolithic variant
# ------------------------------------------------------------

if [ "${RUN_MONOLITHIC}" = true ]; then
    echo "============================================================"
    echo "Cavity flexible bottom FSI regression test (monolithic)"
    echo "End time                    = ${REGRESSION_END_TIME}"
    echo "Final Dy difference        < ${DISP_TOL}"
    echo "Mean Fy difference         < ${FORCE_TOL}"
    echo "============================================================"
    echo

    MONO_DIR="${REGRESSION_ROOT}/monolithic"

    if [ "${CHECK_ONLY}" = false ]; then
        rm -rf "${MONO_DIR}"
        mkdir -p "${MONO_DIR}"
        for item in "${SCRIPT_DIR}"/*; do
            base_item=$(basename "${item}")
            [[ "${base_item}" == "regressionTests" ]] && continue
            cp -a "${item}" "${MONO_DIR}/"
        done

        patch_end_time "${MONO_DIR}/system/controlDict" "${REGRESSION_END_TIME}"

        ( cd "${MONO_DIR}"
          ./Allclean > /dev/null 2>&1 || true
          ./Allrun monolithic > "${ALLRUN_LOGFILE}" 2>&1
        )
    else
        echo "Running in check-only mode: skipping Allclean and Allrun"
    fi

    setup_force_dat "${MONO_DIR}"

    failures=0
    check_results "monolithic" "${MONO_DIR}" || failures=$((failures + $?))

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
fi

# ------------------------------------------------------------
# Partitioned variant
# ------------------------------------------------------------

echo "============================================================"
echo "Cavity flexible bottom FSI regression test"
echo "End time                    = ${REGRESSION_END_TIME}"
echo "Final Dy difference        < ${DISP_TOL}"
echo "Mean Fy difference         < ${FORCE_TOL}"
echo "============================================================"
echo

CASE_DIR="${REGRESSION_ROOT}/main"

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
}

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    patch_end_time "${CASE_DIR}/system/controlDict" "${REGRESSION_END_TIME}"
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

setup_force_dat "${CASE_DIR}"

failures=0
check_results "partitioned" "${CASE_DIR}" || failures=$((failures + $?))

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
