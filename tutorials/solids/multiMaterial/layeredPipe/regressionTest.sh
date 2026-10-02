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

# The point displacement of the removed legacy mechanicalModel, as the max and
# mean component magnitude of its final pointD, from the last commit that had
# it (mcl-stage8-coverage, c3a92b3d). The two are different discretisations -
# the legacy model took its gradient per material on sub-meshes, the framework
# takes it from the material-aware leastSquaresS4f scheme on the whole mesh -
# so their cell displacements differ by about 1.5e-3 of the largest
# displacement, and the point displacements by 1.5e-3 on both foam-extend 4.1
# and OpenFOAM.com v2512. Interpolating the framework displacement to the
# points with foam-extend's least squares fit, which straddles the interface,
# instead puts the difference at 4.9e-3, so the 3e-3 threshold the two were
# compared with separates the two. These are OpenFOAM.com v2512's; foam-extend
# 4.1's differ by 6e-5, well inside it. The legacy model never ran this case
# on OpenFOAM.org, so there it is not compared
REF_POINT_D_MAX=1.62294e-07
REF_POINT_D_MEAN=3.77074821948899e-08
POINT_D_REF_REL_MAX=3e-3

# The parallel arms against the serial arm. The solver is
# not decomposition invariant to round-off with any gradient scheme: a single
# material on this mesh differs from serial by 3.5e-5 in D (interface
# decomposition) and 4.9e-5 with leastSquares, 8.2e-5 with leastSquaresS4f
# (simple), on the legacy path and the framework alike. Two materials measure
# 1.4e-5 on the interface decomposition, where no processor face crosses the
# interface, but 3.5e-4 on the simple one, where the processor boundary cuts
# across it - four times the one-material baseline, on foam-extend 4.1 and
# OpenFOAM.com v2512 alike. That excess is the multi-material path's own
# decomposition dependence where a processor boundary meets the interface,
# tracked in #470; this bound holds it where it is
PARALLEL_REL_MAX=5e-4

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
    # A skip where the application is not built, and a failure in CI, where
    # it always is
    solids4Foam::requireTestApp Test-mechanicalConstitutiveLaw \
        || return $(( $? - 1 ))

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

# The largest pointwise difference between two fields' internal values,
# relative to the first's largest magnitude, and both largest magnitudes, as
# "rel<TAB>max1<TAB>max2"
compare_internal_vector_fields() {
    local a b
    a=$(solids4Foam::internalFieldValues "$1") || return 1
    b=$(solids4Foam::internalFieldValues "$2") || return 1

    paste <(echo "${a}") <(echo "${b}") | awk '
        NF != 2 {
            print "different internalField sizes" > "/dev/stderr"
            bad = 1
            exit 1
        }
        {
            d = $1 - $2; if (d < 0) d = -d; if (d > maxDiff) maxDiff = d
            x = ($1 < 0 ? -$1 : $1); if (x > maxA) maxA = x
            y = ($2 < 0 ? -$2 : $2); if (y > maxB) maxB = y
        }
        END {
            if (bad) exit 1
            printf "%.10g\t%.10g\t%.10g\n", (maxA ? maxDiff/maxA : maxDiff), maxA, maxB
        }'
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
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
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
# Against the legacy answer
# ------------------------------------------------------------
# Two materials sharing an interface, which is the combination the framework
# exists to handle without the legacy per-material subMeshes. The tutorial
# pairs the framework with the material-aware leastSquaresS4f gradient; that
# pairing is the replacement for the subMesh machinery, and routing the
# framework through the subMeshes instead put the radial stress error at
# 0.0305 against the tolerance of 0.03.
#
# The point displacement is where the two part company: the legacy model
# interpolated per material on its sub-meshes, and the framework on the whole
# mesh, extrapolating each cell with its own material's gradient. The stress
# above does not see it, as this solid model does not feed the point
# displacement back into the solution, so it is compared directly
if [ "$CHECK_ONLY" = false ]; then
    if ! grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null
    then
        echo "FAIL: the case constructed no mechanical constitutive law"
        failures=$((failures + 1))
    elif [[ "$(solids4Foam::foamFlavour)" == org ]]; then
        echo "SKIP: pointD against the reference: none on OpenFOAM.org"
    elif ! solids4Foam::checkFieldNorms "pointD" \
        "${CASE_DIR}/$(solids4Foam::latestTime "${CASE_DIR}")/pointD" \
        "${REF_POINT_D_MAX}" "${REF_POINT_D_MEAN}" \
        "${POINT_D_REF_REL_MAX}"
    then
        failures=$((failures + 1))
    fi
fi

# ------------------------------------------------------------
# The framework's multi-material path in parallel
# ------------------------------------------------------------
# Two decompositions: one whose processor boundary is the material interface,
# so an interface face is a processor face on both sides, and a simple one
# whose processor boundaries cut across the interface, so the material
# filter meets processor faces of both materials. Each is held to the
# analytical bound and to the serial arm
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

    ( cd "${dir}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${dir}" && ./Allrun parallel "${decomposition}" \
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
    sr_time="$(solids4Foam::latestTime "${CASE_DIR}")"
    par_time="$(solids4Foam::latestTime "${dir}")"

    if [[ -z "${sr_time}" || "${sr_time}" != "${par_time}" ]]
    then
        echo "FAIL: ${arm}: reached '${par_time}', the serial arm '${sr_time}'"
        failures=$((failures + 1))
        return
    fi

    for fld in D pointD; do
        if ! cmp=$(compare_internal_vector_fields \
            "${CASE_DIR}/${sr_time}/${fld}" "${dir}/${par_time}/${fld}")
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

    ( cd "${dir}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true

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

if [ "$CHECK_ONLY" = false ]; then
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
