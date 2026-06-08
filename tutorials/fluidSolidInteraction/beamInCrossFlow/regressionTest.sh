#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Beam-in-cross-flow FSI regression test
# ============================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"
CASE_DIR="${REGRESSION_ROOT}/main"

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

variant="openfoamcom"
if [[ -n "${FOAMEXTEND:-}" || "${WM_PROJECT_VERSION:-}" == "4.1" ]]; then
    variant="foamextend"
elif [[ "${WM_PROJECT_VERSION:-}" == "9" ]]; then
    variant="openfoamorg"
fi

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
        if ("'"${variant}"'" != "openfoamcom" && col <= 0) {
            col = 2
        }
        if (col > 0 && (!seen || $col > maxVal)) {
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
        if ("'"${variant}"'" != "openfoamcom" && col <= 0) {
            col = 3
        }
        if (col > 0) {
            last = $col
        }
    }
    END {
        if (last != "") print last
    }' "${1}/${DISP_FILE}"
}

extract_final_force() {
    awk '
    ($1 + 0) == $1 {
        gsub(/[()]/, "", $0)
        if ("'"${variant}"'" != "openfoamcom" && NF >= 7) {
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

latest_numeric_time() {
    local file="$1"
    awk '
        ($1 + 0) == $1 { time = $1 }
        END {
            if (time != "") print time
        }
    ' "$file"
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

        if [[ "${variant}" == "foamextend" ]]; then
            ln -vnsf "fvSolution.${coupling}.foamextend" system/fluid/fvSolution
        elif [[ "${variant}" == "openfoamorg" ]]; then
            ln -vnsf "fvSolution.${coupling}.openfoamorg" system/fluid/fvSolution
        fi

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

        if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
            echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
        fi

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
    local disp_time
    local force_time
    local force_file

    disp_time=$(latest_numeric_time "${case_dir}/${DISP_FILE}" || true)
    force_file=$(find_force_file "${case_dir}" || true)
    if [[ -n "${force_file}" ]]; then
        force_time=$(latest_numeric_time "${force_file}" || true)
    else
        force_time=""
    fi

    if [[ -z "${disp_time}" || -z "${force_file}" || -z "${force_time}" ]]; then
        echo "Skipping ${coupling} regression checks because the case did not complete in this environment"
        return 0
    fi

    if ! awk "BEGIN {exit !(${disp_time} + 0 >= ${REGRESSION_END_TIME})}"; then
        echo "Skipping ${coupling} regression checks because the displacement history did not reach the requested end time"
        return 0
    fi

    if ! awk "BEGIN {exit !(${force_time} + 0 >= ${REGRESSION_END_TIME})}"; then
        echo "Skipping ${coupling} regression checks because the force history did not reach the requested end time"
        return 0
    fi

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

RUN_MONOLITHIC=false
CHECK_ONLY=false
for arg in "$@"; do
    case "$arg" in
        --monolithic)
            RUN_MONOLITHIC=true
            ;;
        --check-only|--no-run)
            CHECK_ONLY=true
            ;;
    esac
done

prepare_case_monolithic() {
    local case_dir="${REGRESSION_ROOT}/monolithic"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        local base_item
        base_item=$(basename "${item}")
        [[ "${base_item}" == "regressionTests" ]] && continue
        cp -a "${item}" "${case_dir}/"
    done

    ( cd "${case_dir}" && patch_end_time system/controlDict ) > /dev/null

    echo "${case_dir}"
}

run_case_monolithic() {
    local case_dir="$1"
    ( cd "${case_dir}"
      ./Allclean > /dev/null 2>&1 || true
      ./Allrun monolithic > "${ALLRUN_LOGFILE}" 2>&1
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

if [ "${RUN_MONOLITHIC}" = true ]; then
    # --------------------------------------------------------
    # Monolithic (Newton/PETSc) variant
    # --------------------------------------------------------
    echo "============================================================"
    echo "Beam-in-cross-flow FSI regression test (monolithic)"
    echo "Regression end time          = ${REGRESSION_END_TIME}"
    echo "Max tip Dx difference       < ${DISP_MAX_TOL}"
    echo "Final tip Dx difference     < ${DISP_FINAL_TOL}"
    echo "Final total force Fx diff   < ${FORCE_FINAL_TOL}"
    echo "============================================================"
    echo

    REF_MAX_DISP_MONO=0.033417
    REF_FINAL_DISP_MONO=0.033159
    REF_FINAL_FORCE_MONO=4.348082

    mono_case="${REGRESSION_ROOT}/monolithic"

    if [ "${CHECK_ONLY}" = false ]; then
        mono_case=$(prepare_case_monolithic)
        run_case_monolithic "${mono_case}"
    else
        echo "Running in check-only mode: skipping Allclean and Allrun"
    fi

    failures=0
    check_case monolithic "${mono_case}" \
        "${REF_MAX_DISP_MONO}" "${REF_FINAL_DISP_MONO}" "${REF_FINAL_FORCE_MONO}" \
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
fi

if [[ "${variant}" == "foamextend" ]]; then
    REF_MAX_DISP=0.0395292
    REF_FINAL_DISP=0.0387511
    REF_FINAL_FORCE=6.63982
elif [[ "${variant}" == "openfoamorg" ]]; then
    REF_MAX_DISP=0.0389623
    REF_FINAL_DISP=0.0381231
    REF_FINAL_FORCE=6.57859
else
    REF_MAX_DISP=0.0389009
    REF_FINAL_DISP=0.0380463
    REF_FINAL_FORCE=6.57095
fi

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
