#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# layeredPipe regression test
# Compares sampled transformed stresses against the analytical
# cylinder solution used by the tutorial plots.
# ============================================================

RADIUS_STRESS_ERR_MAX=0.03
THETA_POINT_ERR_MAX=0.01

# Largest relative difference allowed between the point displacements of the
# legacy and framework arms. The two are different discretisations - the legacy
# arm takes its gradient per material on sub-meshes, the framework arm from the
# material-aware leastSquaresS4f scheme on the whole mesh - so their cell
# displacements differ by about 1.5e-3 of the largest displacement, and the
# point displacements by 1.5e-3 on both foam-extend 4.1 and OpenFOAM.com
# v2512. Interpolating the framework displacement to the points with
# foam-extend's least squares fit, which straddles the interface, instead
# puts the difference at 4.9e-3, so the threshold separates the two
POINT_D_ARM_REL_MAX=3e-3

# The parallel framework arms against the serial framework arm. The solver is
# not decomposition invariant to round-off with any gradient scheme: a single
# material on this mesh, decomposed as the interface arm is, differs from
# serial by 3.5e-5 in D on the legacy path with leastSquares and on the
# framework with leastSquaresS4f alike, identically on foam-extend 4.1 and
# OpenFOAM.com v2512. The two materials measure 1.4e-5 (interface) and 9e-5
# (simple), the same order, so this bounds the decomposition dependence of
# the multi-material machinery rather than asserting there is none
PARALLEL_REL_MAX=2e-4

R1=0.05
R2=0.07
R3=0.1
E1=20e9
E2=200e9
NU1=0.35
NU2=0.3
PLOAD=1e5

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
CONSTITUTIVE_LOGFILE="log.Test-mechanicalConstitutiveLaw"

echo "============================================================"
echo "layeredPipe regression test"
echo "Max radial stress error  < ${RADIUS_STRESS_ERR_MAX}"
echo "Outer-point theta error  < ${THETA_POINT_ERR_MAX}"
echo "Plus the mechanicalConstitutiveLaw framework checks"
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
}

sample_file() {
    local preferred
    preferred="$(find "${CASE_DIR}" -name 'line_sigma:Transformed.xy' | sort | tail -n 1)"
    if [[ -n "${preferred}" ]]; then
        printf '%s\n' "${preferred}"
        return 0
    fi

    find "${CASE_DIR}" -name 'sigma:Transformed' | sort | tail -n 1
}

# Exercise the mechanicalConstitutiveLawManager evaluation paths on this case.
# This tutorial is used because it is the only one with more than one material,
# and the manager's integration-point addressing is per material. The solid
# models can now take their stress from the framework, so this is no longer its
# only runtime coverage, but it remains the only multi-material coverage
run_constitutive_test() {
    if ! command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
        echo "SKIP: Test-mechanicalConstitutiveLaw not found in PATH"
        return 0
    fi

    if [[ ! -d "${CASE_DIR}/constant/polyMesh" ]]; then
        echo "SKIP: mechanicalConstitutiveLaw checks (case has no mesh)"
        return 0
    fi

    if ( cd "${CASE_DIR}" && Test-mechanicalConstitutiveLaw \
            > "${CONSTITUTIVE_LOGFILE}" 2>&1 )
    then
        local n_passed
        n_passed=$(grep -c 'PASS:' "${CASE_DIR}/${CONSTITUTIVE_LOGFILE}" || true)

        if (( n_passed == 0 )); then
            # The test exited zero but reported nothing: it did not run the
            # checks it is here for, so this is a failure, not a skip
            echo "FAIL: mechanicalConstitutiveLaw checks (no checks reported)"
            return 1
        fi

        echo "PASS: mechanicalConstitutiveLaw checks (${n_passed} checks)"
        return 0
    fi

    echo "FAIL: mechanicalConstitutiveLaw checks"
    grep 'FAIL:' "${CASE_DIR}/${CONSTITUTIVE_LOGFILE}" || true
    return 1
}

# Relative difference between the internal fields of two vector fields on the
# same mesh, followed by the largest magnitude component of each, separated by
# tabs as IFS does not split on spaces here
compare_internal_vector_fields() {
    python3 - "$1" "$2" << 'PYEOF'
import re
import sys

number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

def read_internal(path):
    text = open(path).read()
    nonuniform = re.search(
        r"\binternalField\s+nonuniform\s+List<vector>\s+\d+\s*\((.*?)\)\s*;",
        text,
        re.DOTALL,
    )
    if not nonuniform:
        raise ValueError(f"cannot parse internalField in {path}")

    values = re.findall(
        rf"\(({number})\s+({number})\s+({number})\)", nonuniform.group(1)
    )
    if not values:
        raise ValueError(f"empty internalField in {path}")
    return [tuple(map(float, value)) for value in values]

try:
    a = read_internal(sys.argv[1])
    b = read_internal(sys.argv[2])
    if len(a) != len(b):
        raise ValueError("different internalField sizes")
    max_diff = max(abs(x - y) for av, bv in zip(a, b) for x, y in zip(av, bv))
    max_a = max(abs(x) for av in a for x in av)
    max_b = max(abs(x) for av in b for x in av)
    print(f"{max_diff/max_a if max_a else max_diff:.10g}\t{max_a:.10g}\t{max_b:.10g}")
except (OSError, ValueError) as error:
    print(error, file=sys.stderr)
    sys.exit(1)
PYEOF
}

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

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun "${ARM:-}" > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

DATA_FILE="$(sample_file)"

if [[ -z "${DATA_FILE}" || ! -f "${DATA_FILE}" ]]; then
    echo "FAIL: Could not find sampled stress data"
    exit 1
fi

compute_radial_err() {
    awk -v r1="${R1}" -v r2="${R2}" -v r3="${R3}" \
        -v e1="${E1}" -v e2="${E2}" -v nu1="${NU1}" -v nu2="${NU2}" \
        -v p="${PLOAD}" '
    function abs(x) { return x < 0 ? -x : x }
    BEGIN {
        pint = (2*r1*r1*p/(e1*(r2*r2-r1*r1))) / (((1.0/e2)*(((r3*r3+r2*r2)/(r3*r3-r2*r2))+nu2)) + ((1.0/e1)*(((r2*r2+r1*r1)/(r2*r2-r1*r1))-nu1)))
        maxRad = 0
    }
    NF >= 5 {
        r = $1
        sigmaR = $2

        if (r < r2) {
            sigmaRAnal = (r1*r1*p - r2*r2*pint + (pint-p)*(r1*r2/r)^2) / (r2*r2-r1*r1)
        } else {
            sigmaRAnal = (r2*r2*pint - pint*(r2*r3/r)^2) / (r3*r3-r2*r2)
        }

        radialErr = abs(sigmaR - sigmaRAnal) / p
        if (radialErr > maxRad) {
            maxRad = radialErr
        }
    }
    END {
        printf "%.12g\n", maxRad
    }
    ' "$1"
}

max_radial_err="$(compute_radial_err "${DATA_FILE}")"

outer_theta_err="$(awk -v r1="${R1}" -v r2="${R2}" -v r3="${R3}" \
        -v e1="${E1}" -v e2="${E2}" -v nu1="${NU1}" -v nu2="${NU2}" -v p="${PLOAD}" '
    function abs(x) { return x < 0 ? -x : x }
    BEGIN {
        pint = (2*r1*r1*p/(e1*(r2*r2-r1*r1))) / (((1.0/e2)*(((r3*r3+r2*r2)/(r3*r3-r2*r2))+nu2)) + ((1.0/e1)*(((r2*r2+r1*r1)/(r2*r2-r1*r1))-nu1)))
        lastThetaAnal = 0
        lastThetaSample = 0
    }
    NF >= 5 {
        r = $1
        sigmaTheta = $5

        if (r < r2) {
            sigmaThetaAnal = (r1*r1*p - r2*r2*pint - (pint-p)*(r1*r2/r)^2) / (r2*r2-r1*r1)
        } else {
            sigmaThetaAnal = (r2*r2*pint + pint*(r2*r3/r)^2) / (r3*r3-r2*r2)
        }

        lastThetaAnal = sigmaThetaAnal
        lastThetaSample = sigmaTheta
    }
    END {
        thetaErr = abs(lastThetaSample - lastThetaAnal) / p
        printf "%.12g\n", thetaErr
    }
    ' "${DATA_FILE}")"

if [[ -z "${max_radial_err}" || -z "${outer_theta_err}" ]]; then
    echo "FAIL: Could not extract regression errors"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${max_radial_err} < ${RADIUS_STRESS_ERR_MAX})}"; then
    printf "PASS: Max radial stress error = %.6g\n" "${max_radial_err}"
else
    printf "FAIL: Max radial stress error = %.6g\n" "${max_radial_err}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${outer_theta_err} < ${THETA_POINT_ERR_MAX})}"; then
    printf "PASS: Outer-point theta error = %.6g\n" "${outer_theta_err}"
else
    printf "FAIL: Outer-point theta error = %.6g\n" "${outer_theta_err}"
    failures=$((failures + 1))
fi

if ! run_constitutive_test; then
    failures=$((failures + 1))
fi

# ------------------------------------------------------------
# The same case on the constitutive-law framework
# ------------------------------------------------------------
# Two materials sharing an interface, which is the combination the framework
# exists to handle without the legacy per-material subMeshes. The framework
# arm pairs the switch with the material-aware leastSquaresS4f gradient; that
# pairing is the replacement for the subMesh machinery, and routing the
# framework through the subMeshes instead puts the radial stress error at
# 0.0305 against a tolerance of 0.03
if [ "$CHECK_ONLY" = false ]; then
    FRAMEWORK_DIR="${REGRESSION_ROOT}/framework"
    rm -rf "${FRAMEWORK_DIR}"; mkdir -p "${FRAMEWORK_DIR}"
    for item in "${SCRIPT_DIR}"/*; do
        [[ "$(basename "${item}")" == "regressionTests" ]] && continue
        cp -a "${item}" "${FRAMEWORK_DIR}/"
    done

    ( cd "${FRAMEWORK_DIR}" && ./Allrun framework > "${ALLRUN_LOGFILE}" 2>&1 ) \
        || true

    if solids4Foam::regressionCaseSkipped "${FRAMEWORK_DIR}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: the framework arm does not run in this environment"
    elif ! grep -q "Selecting mechanical constitutive law" \
        "${FRAMEWORK_DIR}/log.solids4Foam" 2>/dev/null
    then
        echo "FAIL: the framework arm did not use the framework"
        failures=$((failures + 1))
    else
        fw_file="$(find "${FRAMEWORK_DIR}" -name 'line_sigma:Transformed.xy' \
            | sort | tail -n 1)"

        if [[ -z "${fw_file}" ]]; then
            echo "FAIL: the framework arm produced no sampled stress"
            failures=$((failures + 1))
        else
            fw_radial="$(compute_radial_err "${fw_file}")"

            if awk "BEGIN {exit !(${fw_radial} < ${RADIUS_STRESS_ERR_MAX})}"
            then
                printf "PASS: framework: Max radial stress error = %.6g\n" \
                    "${fw_radial}"
            else
                printf "FAIL: framework: Max radial stress error = %.6g\n" \
                    "${fw_radial}"
                failures=$((failures + 1))
            fi
        fi

        # The point displacement is where the two arms part company: the
        # legacy arm interpolates per material on its sub-meshes, and the
        # framework arm on the whole mesh, extrapolating each cell with its
        # own material's gradient. The stress above does not see it, as this
        # solid model does not feed the point displacement back into the
        # solution, so compare it directly against the legacy arm
        lg_time="$(solids4Foam::latestTime "${CASE_DIR}")"
        fw_time="$(solids4Foam::latestTime "${FRAMEWORK_DIR}")"
        lg_pointD="${CASE_DIR}/${lg_time}/pointD"
        fw_pointD="${FRAMEWORK_DIR}/${fw_time}/pointD"

        if grep -q "Selecting mechanical constitutive law" \
            "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null
        then
            # Only when ARM selects the framework for the main case too
            echo "SKIP: point displacement comparison (main arm is not legacy)"
        elif [[ -z "${lg_time}" || "${lg_time}" != "${fw_time}" ]]; then
            echo "FAIL: the arms reached different times" \
                "('${lg_time}' vs '${fw_time}')"
            failures=$((failures + 1))
        elif [[ ! -f "${lg_pointD}" || ! -f "${fw_pointD}" ]]; then
            echo "FAIL: an arm wrote no pointD"
            failures=$((failures + 1))
        elif ! pointD_cmp=$(compare_internal_vector_fields \
            "${lg_pointD}" "${fw_pointD}")
        then
            echo "FAIL: could not compare the pointD fields"
            failures=$((failures + 1))
        else
            read -r pointD_rel lg_max fw_max <<< "${pointD_cmp}"

            # A zero legacy field would make any comparison pass
            if ! awk "BEGIN {exit !(${lg_max} > 1e-9 && ${fw_max} > 1e-9)}"
            then
                printf "FAIL: pointD is trivially small (%.4g, %.4g)\n" \
                    "${lg_max}" "${fw_max}"
                failures=$((failures + 1))
            elif awk "BEGIN {exit !(${pointD_rel} < ${POINT_D_ARM_REL_MAX})}"
            then
                printf "PASS: framework: pointD relative diff to legacy = %.4g\n" \
                    "${pointD_rel}"
            else
                printf "FAIL: framework: pointD relative diff to legacy = %.4g\n" \
                    "${pointD_rel}"
                failures=$((failures + 1))
            fi
        fi
    fi
fi

# ------------------------------------------------------------
# The framework's multi-material path in parallel
# ------------------------------------------------------------
# Two decompositions: one whose processor boundary is the material interface,
# so an interface face is a processor face on both sides, and a simple one
# whose processor boundaries cut across the interface, so the material
# filter meets processor faces of both materials. Each is held to the
# analytical bound and to the serial framework arm
run_parallel_arm() {
    local decomposition="$1"
    local arm="parallel-${decomposition}"
    local dir="${REGRESSION_ROOT}/${arm}"
    local item

    rm -rf "${dir}"; mkdir -p "${dir}"
    for item in "${SCRIPT_DIR}"/*; do
        [[ "$(basename "${item}")" == "regressionTests" ]] && continue
        cp -a "${item}" "${dir}/"
    done

    ( cd "${dir}" && ./Allrun framework parallel "${decomposition}" \
        > "${ALLRUN_LOGFILE}" 2>&1 ) || true

    if solids4Foam::regressionCaseSkipped "${dir}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: ${arm} does not run in this environment"
        return
    fi

    if ! grep -q "Selecting mechanical constitutive law" \
        "${dir}/${SOLVER_LOGFILE}" 2>/dev/null
    then
        echo "FAIL: ${arm}: did not use the framework"
        failures=$((failures + 1))
        return
    fi

    if [[ ! -d "${dir}/processor1" ]]
    then
        echo "FAIL: ${arm}: did not run in parallel"
        failures=$((failures + 1))
        return
    fi

    # The interface decomposition is only the test it claims to be if
    # processor 0 holds exactly the inner material
    if [[ "${decomposition}" == "interface" ]]
    then
        local nInner nProc0
        nInner=$(sed -n '/^[0-9][0-9]*$/{p;q}' \
            "${dir}/constant/polyMesh/sets/inner")
        nProc0=$(grep -o 'nCells:[ ]*[0-9]*' \
            "${dir}/processor0/constant/polyMesh/owner" | grep -o '[0-9]*$')
        if [[ -z "${nInner}" || "${nInner}" != "${nProc0}" ]]
        then
            echo "FAIL: ${arm}: processor 0 holds ${nProc0} cells," \
                "the inner material ${nInner}"
            failures=$((failures + 1))
            return
        fi
    fi

    local par_file par_radial
    par_file="$(find "${dir}" -name 'line_sigma:Transformed.xy' \
        -not -path '*/processor*' | sort | tail -n 1)"

    if [[ -z "${par_file}" ]]
    then
        echo "FAIL: ${arm}: produced no sampled stress"
        failures=$((failures + 1))
    else
        par_radial="$(compute_radial_err "${par_file}")"
        if awk "BEGIN {exit !(${par_radial} < ${RADIUS_STRESS_ERR_MAX})}"
        then
            printf "PASS: %s: Max radial stress error = %.6g\n" \
                "${arm}" "${par_radial}"
        else
            printf "FAIL: %s: Max radial stress error = %.6g\n" \
                "${arm}" "${par_radial}"
            failures=$((failures + 1))
        fi
    fi

    local sr_time par_time fld cmp rel sr_max par_max
    sr_time="$(solids4Foam::latestTime "${FRAMEWORK_DIR}")"
    par_time="$(solids4Foam::latestTime "${dir}")"

    if [[ -z "${sr_time}" || "${sr_time}" != "${par_time}" ]]
    then
        echo "FAIL: ${arm}: reached '${par_time}', the serial arm '${sr_time}'"
        failures=$((failures + 1))
        return
    fi

    for fld in D pointD; do
        if ! cmp=$(compare_internal_vector_fields \
            "${FRAMEWORK_DIR}/${sr_time}/${fld}" "${dir}/${par_time}/${fld}")
        then
            echo "FAIL: ${arm}: could not compare ${fld} with the serial arm"
            failures=$((failures + 1))
            continue
        fi

        read -r rel sr_max par_max <<< "${cmp}"

        if ! awk "BEGIN {exit !(${sr_max} > 1e-9 && ${par_max} > 1e-9)}"
        then
            printf "FAIL: %s: %s is trivially small (%.4g, %.4g)\n" \
                "${arm}" "${fld}" "${sr_max}" "${par_max}"
            failures=$((failures + 1))
        elif awk "BEGIN {exit !(${rel} < ${PARALLEL_REL_MAX})}"
        then
            printf "PASS: %s: %s relative diff to serial = %.4g\n" \
                "${arm}" "${fld}" "${rel}"
        else
            printf "FAIL: %s: %s relative diff to serial = %.4g\n" \
                "${arm}" "${fld}" "${rel}"
            failures=$((failures + 1))
        fi
    done
}

# ------------------------------------------------------------
# A material one cell thick is refused, and for that reason
# ------------------------------------------------------------
# The inner material is reduced to the single layer of cells between r = 70
# and 71 mm, inside the outer one. Its cells then have no neighbour of their
# own material in the radial direction, even after the stencil is widened to
# point neighbours, so there is no gradient to reconstruct and the
# material-aware scheme must stop, naming the rank-deficient cells, rather
# than build one from both materials. The error's other cause, a stencil
# truncated at a processor boundary, is not reachable on this mesh: the
# compact stencil reaches across processor faces, and a processor slab one
# layer thick runs to completion
run_one_cell_thick_arm() {
    local arm="oneCellThick"
    local dir="${REGRESSION_ROOT}/${arm}"
    local item

    rm -rf "${dir}"; mkdir -p "${dir}"
    for item in "${SCRIPT_DIR}"/*; do
        [[ "$(basename "${item}")" == "regressionTests" ]] && continue
        cp -a "${item}" "${dir}/"
    done

    cat > "${dir}/batch.setSet" << 'SETEOF'
cellSet outer new cylinderToCell (0.0 0.0 -100) (0.0 0.0 100) 100e-3
cellSet inner new cylinderToCell (0.0 0.0 -100) (0.0 0.0 100) 71e-3
cellSet core new cylinderToCell (0.0 0.0 -100) (0.0 0.0 100) 70e-3
cellSet inner delete cellToCell core
cellSet outer delete cellToCell inner
SETEOF

    ( cd "${dir}" && ./Allrun framework > "${ALLRUN_LOGFILE}" 2>&1 ) || true

    if solids4Foam::regressionCaseSkipped "${dir}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: ${arm} does not run in this environment"
        return
    fi

    local nLayer
    nLayer=$(sed -n '/^[0-9][0-9]*$/{p;q}' \
        "${dir}/constant/polyMesh/sets/inner" 2>/dev/null)

    if grep -q "^End" "${dir}/${SOLVER_LOGFILE}" 2>/dev/null
    then
        echo "FAIL: ${arm}: ran to completion on a material one cell thick"
        failures=$((failures + 1))
    elif grep -q "${nLayer} cells remain rank-deficient after widening the gradient stencil to point neighbours of the same material" \
        "${dir}/${SOLVER_LOGFILE}" 2>/dev/null
    then
        echo "PASS: ${arm}: refused, naming the ${nLayer} cells of the layer"
    else
        echo "FAIL: ${arm}: did not stop for the rank-deficient stencil" \
            "of the ${nLayer} layer cells"
        failures=$((failures + 1))
    fi
}

if [ "$CHECK_ONLY" = false ] && [[ -n "${FRAMEWORK_DIR:-}" ]] \
    && ! solids4Foam::regressionCaseSkipped "${FRAMEWORK_DIR}/${ALLRUN_LOGFILE}"
then
    run_parallel_arm interface
    run_parallel_arm simple
    run_one_cell_thick_arm
fi

if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
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
