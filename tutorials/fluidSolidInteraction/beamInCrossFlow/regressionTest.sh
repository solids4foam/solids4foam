#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Beam-in-cross-flow FSI regression test
# ============================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

# Relative tolerances. All supported OpenFOAM versions share one set of
# references (the high-order solid has its own), so these cover the small
# differences between coupling methods (Aitken, IQN-ILS, Robin: < 3e-4),
# OpenFOAM versions (< 1e-4), compilers and platforms.
DISP_MAX_TOL=1e-2
DISP_FINAL_TOL=1e-2
FORCE_FINAL_TOL=5e-3

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_displacement.dat"

# Shortened end time for regression efficiency
REGRESSION_END_TIME=1

# The diffusionElastic fluid mesh motion variant costs about ten times as much
# as the others, so it is run to a shorter end time, with its own references
# and tolerances
SOLID_MESH_MOTION_END_TIME=0.5
SOLID_MESH_MOTION_DISP_TOL=1e-3
SOLID_MESH_MOTION_FORCE_TOL=0.1

# foam-extend is not supported for this case: the Allrun declares a skip
variant="openfoamcom"
if [[ "${WM_PROJECT_VERSION:-}" != v* ]]; then
    variant="openfoamorg"
fi

echo "============================================================"
echo "Beam-in-cross-flow FSI regression test"
echo "Regression end time          = ${REGRESSION_END_TIME}"
echo "Max point-A Dx rel. diff.    < ${DISP_MAX_TOL}"
echo "Final point-A Dx rel. diff.  < ${DISP_FINAL_TOL}"
echo "Final total Fx rel. diff.    < ${FORCE_FINAL_TOL}"
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
    }' "${1}"
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

# |value - reference| / |reference|
rel_diff() {
    awk -v x="$1" -v r="$2" 'BEGIN {d = x - r; if (d < 0) d = -d; if (r < 0) r = -r; print d/r}'
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
    local case_name="${2:-${coupling}}"
    local case_dir="${REGRESSION_ROOT}/${case_name}"
    local solution_setup="${coupling}"
    local interface_condition="dirichletNeumann"

    if [[ "${coupling}" == "robin" ]]; then
        solution_setup="iqnils"
        interface_condition="robin"
    fi

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        local base_item
        base_item=$(basename "${item}")

        if [[ "${base_item}" == "regressionTests" || "${base_item}" == "verification" ]]; then
            continue
        fi

        cp -a "${item}" "${case_dir}/"
    done

    (
        cd "${case_dir}"

        ln -vnsf "fsiProperties.${coupling}" constant/fsiProperties
        ln -vnsf "solidProperties.${solution_setup}" constant/solid/solidProperties
        ln -vnsf "controlDict.${solution_setup}" system/controlDict
        ln -vnsf "fvSolution.${solution_setup}" system/fluid/fvSolution
        ln -vnsf "U.${interface_condition}" 0/fluid/U
        ln -vnsf "p.${interface_condition}" 0/fluid/p

        if [[ "${variant}" == "openfoamorg" ]]; then
            ln -vnsf "fvSolution.${solution_setup}.openfoamorg" system/fluid/fvSolution
        fi

        patch_end_time "system/controlDict.${solution_setup}"
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

run_mesh_motion_case() {
    local case_dir="$1"
    local mesh_motion="$2"

    echo "Running IQNILS ${mesh_motion} fluid mesh motion regression case"
    (
        cd "${case_dir}"
        ./Allclean > /dev/null 2>&1 || true
        ./Allrun iqnils "${mesh_motion}" > "${ALLRUN_LOGFILE}" 2>&1
    )
}

run_high_order_case() {
    local case_dir="$1"

    echo "Running high-order IQNILS regression case"
    (
        cd "${case_dir}"
        ./Allclean > /dev/null 2>&1 || true
        ./Allrun iqnils highOrder > "${ALLRUN_LOGFILE}" 2>&1
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
    local disp_diff_abs
    local final_disp_diff_abs
    local final_force_diff_abs
    local failures=0
    local disp_time
    local force_time
    local force_file

    # A skip is only valid if the tutorial declared one in the Allrun log.
    # Anything else that leaves the expected output missing or incomplete is a
    # failure.
    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "Skipping ${coupling} regression checks because the tutorial skipped in this environment"
        return 0
    fi

    disp_time=$(latest_numeric_time "${case_dir}/${DISP_FILE}" || true)
    force_file=$(find_force_file "${case_dir}" || true)
    if [[ -n "${force_file}" ]]; then
        force_time=$(latest_numeric_time "${force_file}" || true)
    else
        force_time=""
    fi

    if [[ -z "${disp_time}" || -z "${force_file}" || -z "${force_time}" ]]; then
        echo "FAIL [${coupling}]: the case did not run or did not complete in this"
        echo "      environment: expected output is missing and the tutorial did"
        echo "      not declare a skip (see ${case_dir}/${ALLRUN_LOGFILE})"
        return 1
    fi

    if ! awk "BEGIN {exit !(${disp_time} + 0 >= ${REGRESSION_END_TIME})}"; then
        echo "FAIL [${coupling}]: the displacement history stops at t = ${disp_time},"
        echo "      short of the requested end time ${REGRESSION_END_TIME}: the case"
        echo "      did not complete"
        return 1
    fi

    if ! awk "BEGIN {exit !(${force_time} + 0 >= ${REGRESSION_END_TIME})}"; then
        echo "FAIL [${coupling}]: the force history stops at t = ${force_time},"
        echo "      short of the requested end time ${REGRESSION_END_TIME}: the case"
        echo "      did not complete"
        return 1
    fi

    max_disp=$(extract_max_displacement "${case_dir}")
    final_disp=$(extract_final_displacement "${case_dir}")
    final_force=$(extract_final_force "${force_file}")

    if [[ -z "${max_disp}" || -z "${final_disp}" || -z "${final_force}" ]]; then
        echo "FAIL [${coupling}]: Could not extract regression quantities"
        return 1
    fi

    disp_diff_abs=$(rel_diff "${max_disp}" "${ref_max_disp}")

    final_disp_diff_abs=$(rel_diff "${final_disp}" "${ref_final_disp}")

    final_force_diff_abs=$(rel_diff "${final_force}" "${ref_final_force}")

    if awk "BEGIN {exit !(${disp_diff_abs} < ${DISP_MAX_TOL})}"; then
        printf "PASS [%s]: max displacement = %.6g (rel. Δ = %.3g)\n" \
            "${coupling}" "${max_disp}" "${disp_diff_abs}"
    else
        printf "FAIL [%s]: max displacement = %.6g (rel. Δ = %.3g)\n" \
            "${coupling}" "${max_disp}" "${disp_diff_abs}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${final_disp_diff_abs} < ${DISP_FINAL_TOL})}"; then
        printf "PASS [%s]: final tip Dx = %.6g (rel. Δ = %.3g)\n" \
            "${coupling}" "${final_disp}" "${final_disp_diff_abs}"
    else
        printf "FAIL [%s]: final tip Dx = %.6g (rel. Δ = %.3g)\n" \
            "${coupling}" "${final_disp}" "${final_disp_diff_abs}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${final_force_diff_abs} < ${FORCE_FINAL_TOL})}"; then
        printf "PASS [%s]: final total Fx = %.6g (rel. Δ = %.3g)\n" \
            "${coupling}" "${final_force}" "${final_force_diff_abs}"
    else
        printf "FAIL [%s]: final total Fx = %.6g (rel. Δ = %.3g)\n" \
            "${coupling}" "${final_force}" "${final_force_diff_abs}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

# Reference values at t = 1 s for the tutorial's Richter definition (inlet
# peak 0.2 m/s, St Venant-Kirchhoff E = 1.4 MPa), from the IQN-ILS runs of the
# GitHub CI OpenFOAM-v2512-PETSc image. One set serves all supported versions:
# the OpenFOAM-9-PETSc image agrees to within 4e-5. The high-order solid gives
# a different, less stiff answer on this coarse mesh and therefore has its own
# references.
REF_MAX_DISP=8.27140e-5
REF_FINAL_DISP=4.41995e-5
REF_FINAL_FORCE=1.40949
REF_HO_MAX_DISP=9.64144e-5
REF_HO_FINAL_DISP=5.22025e-5
REF_HO_FINAL_FORCE=1.41024

# diffusionElastic fluid mesh motion, to t = 0.5, from the same CI image (it
# runs only on OpenFOAM.com): the tip displacement is still increasing, so its
# maximum is its final value
REF_SOLID_MESH_MOTION_MAX_DISP=7.76066e-5
REF_SOLID_MESH_MOTION_FINAL_DISP=7.76066e-5
REF_SOLID_MESH_MOTION_FINAL_FORCE=3.00683

aitken_case=$(prepare_case aitken)
run_case "${aitken_case}" aitken

iqnils_case=$(prepare_case iqnils)
run_case "${iqnils_case}" iqnils

robin_case=$(prepare_case robin)
run_case "${robin_case}" robin

high_order_case=$(prepare_case iqnils highOrder)
run_high_order_case "${high_order_case}"

# The fluid mesh motion variants: displacementLaplacian changes the beam
# response by much less than the tolerances, so it shares the references;
# diffusionElastic is run to the shorter end time, and diffusionHyperElastic
# is not run. Each Allrun declares a skip where its mesh motion is unavailable
displacement_case=$(prepare_case iqnils displacementLaplacian)
run_mesh_motion_case "${displacement_case}" displacementLaplacian

solid_mesh_motion_case=$(
    REGRESSION_END_TIME="${SOLID_MESH_MOTION_END_TIME}" \
        prepare_case iqnils diffusionElastic
)
run_mesh_motion_case "${solid_mesh_motion_case}" diffusionElastic

echo

failures=0

check_case aitken "${aitken_case}" \
    "${REF_MAX_DISP}" "${REF_FINAL_DISP}" "${REF_FINAL_FORCE}" \
    || failures=$((failures + $?))

check_case iqnils "${iqnils_case}" \
    "${REF_MAX_DISP}" "${REF_FINAL_DISP}" "${REF_FINAL_FORCE}" \
    || failures=$((failures + $?))

check_case robin "${robin_case}" \
    "${REF_MAX_DISP}" "${REF_FINAL_DISP}" "${REF_FINAL_FORCE}" \
    || failures=$((failures + $?))

check_case highOrder "${high_order_case}" \
    "${REF_HO_MAX_DISP}" "${REF_HO_FINAL_DISP}" "${REF_HO_FINAL_FORCE}" \
    || failures=$((failures + $?))

check_case displacementLaplacian "${displacement_case}" \
    "${REF_MAX_DISP}" "${REF_FINAL_DISP}" "${REF_FINAL_FORCE}" \
    || failures=$((failures + $?))

REGRESSION_END_TIME="${SOLID_MESH_MOTION_END_TIME}" \
DISP_MAX_TOL="${SOLID_MESH_MOTION_DISP_TOL}" \
DISP_FINAL_TOL="${SOLID_MESH_MOTION_DISP_TOL}" \
FORCE_FINAL_TOL="${SOLID_MESH_MOTION_FORCE_TOL}" \
    check_case diffusionElastic "${solid_mesh_motion_case}" \
    "${REF_SOLID_MESH_MOTION_MAX_DISP}" "${REF_SOLID_MESH_MOTION_FINAL_DISP}" \
    "${REF_SOLID_MESH_MOTION_FINAL_FORCE}" \
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
