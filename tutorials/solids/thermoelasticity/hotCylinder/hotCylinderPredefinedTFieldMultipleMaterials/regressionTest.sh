#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# GNU sed, for the in-place edits below
solids4Foam::requireGnuSed

# ============================================================
# hotCylinderPredefinedTFieldMultipleMaterials regression test
#
# Two materials, steel and aluminium, whose temperature is not solved for but
# read at each time step from the separate case hotCylinderTemperatureField
# (TcaseDirectory). Two arms, each held to the answer of the removed legacy
# mechanicalModel:
#
#   main              the tutorial as it is
#   single            steel alone
#
# The single-material arm isolates the reading of T. With one material the
# framework discretised the problem as the legacy model did, so any difference
# beyond the solution tolerance would be the temperature it was given. The
# two-material arm is then the multi-material path with that input.
# ============================================================

# Bounds on the main arm's final extrema. The temperature is uniform in
# each time directory, so the stress is the mismatch of the two materials'
# expansion. Measured on the legacy model: max epsilonEq 1.47489e-3 and max
# sigmaEq 1.8164e8 on OpenFOAM.com v2512 and OpenFOAM.org 9, 1.47582e-3 and
# 1.81641e8 on foam-extend 4.1
EPS_MIN=1.4e-3
EPS_MAX=1.55e-3
SIGMA_MIN=1.75e8
SIGMA_MAX=1.9e8

# The final D of the removed legacy mechanicalModel, as the max and mean
# component magnitude of the field written to fourteen figures, for each arm,
# from the last commit that had it (mcl-stage8-coverage, c3a92b3d). These are
# OpenFOAM.com v2512's; OpenFOAM.org 9 agrees to 6e-8, and foam-extend 4.1
# differs by 1.7e-4 in the two-material arm and 2e-5 in the single-material one
REF_D_MAX=0.00098799181372038
REF_D_MEAN=0.000340730093870825
REF_SINGLE_D_MAX=0.0007175233920732
REF_SINGLE_D_MEAN=0.000264826849529889

# Two materials. These are different discretisations of the interface - per
# material sub-meshes on the legacy path, the material-aware leastSquaresS4f
# gradient on one mesh on the framework - as on layeredPipe and punch, so they
# agree only to the discretisation. Measured: 1.23e-3 of the largest
# displacement on OpenFOAM.com v2512 and OpenFOAM.org 9, 1.35e-3 on
# foam-extend 4.1. This is the coarse check: a 1% change in steel's alpha
# took the framework to 6.2e-3 and failed, a 0.1% change (1.7e-3) did not, and
# holding T at its second-step value took it to 0.12. The single-material arm
# below is the fine one
FRAMEWORK_D_REL_TOL=3e-3

# One material. The same discretisation, so the two agreed to the solution
# tolerance (1e-6): measured 2.7e-8 on OpenFOAM.com v2512, 3.1e-8 on
# foam-extend 4.1 and 8.2e-8 on OpenFOAM.org 9. The threshold, 5e-5, covers
# the 2e-5 between the forks and is below the 1e-4 that a 0.01% change in
# alpha makes
SINGLE_D_REL_TOL=5e-5

# The last time step, and the last directory of the temperature case
COMPARISON_END_TIME=4

# A thermal stress of at least this much must be present at the end in every
# arm, so that the comparisons are not of stress-free solutions
SIGMA_NONZERO_MIN=1e7

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
FRAMEWORK_MARK="Selecting mechanical constitutive law"
FRAMEWORK_T_READ="Reading T for the mechanical constitutive laws from"

echo "============================================================"
echo "hotCylinderPredefinedTFieldMultipleMaterials regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "D against the reference, relative to its largest value: two"\
" materials < ${FRAMEWORK_D_REL_TOL}, one material < ${SINGLE_D_REL_TOL}"
echo "============================================================"
echo

copy_case() {
    local dest="$1"

    rm -rf "${dest}"
    mkdir -p "${dest}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${dest}/"
    done

    # Enough digits that the comparisons measure the solutions rather than
    # the last digit written
    "${SOLIDS4FOAM_SED}" -i 's|^writePrecision.*|writePrecision  14;|' \
        "${dest}/system/controlDict"
}

# Keep steel only. One law covers the whole mesh, so no cellZone is used
use_single_material() {
    local dir="$1"

    # From the aluminium entry up to the list's closing bracket
    "${SOLIDS4FOAM_SED}" -i '/^    aluminium/,/^);/{/^);/!d}' \
        "${dir}/constant/mechanicalProperties"

    if grep -q aluminium "${dir}/constant/mechanicalProperties"; then
        echo "FAIL: could not remove the second material from ${dir}"
        return 1
    fi

    # One material needs no material-aware gradient, and the legacy single
    # material answer this arm is held to was computed with the gradient the
    # tutorial used before it needed one. Using it here keeps the comparison
    # about the temperature the law was given, not about the gradient
    "${SOLIDS4FOAM_SED}" -i \
        's|^\( *default *\)leastSquaresS4f;|\1pointCellsLeastSquares;|' \
        "${dir}/system/fvSchemes"

    if ! grep -q "pointCellsLeastSquares" "${dir}/system/fvSchemes"; then
        echo "FAIL: could not set the single-material gradient in ${dir}"
        return 1
    fi
}

extract_last() {
    grep "$1" "$2/${SOLVER_LOGFILE}" | tail -n 1 | awk '{print $NF}'
}

# Checks common to every arm: it completed, took its material and its
# temperature from the framework, and ended at the last time with a real
# thermal stress
check_arm() {
    local dir="$1"

    if ! grep -q "^End" "${dir}/${SOLVER_LOGFILE}" 2>/dev/null \
      || grep -qE "FOAM FATAL|did not converge" "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${dir##*/} did not complete and converge"
        return 1
    fi

    if ! grep -q "${FRAMEWORK_MARK}" "${dir}/${SOLVER_LOGFILE}"; then
        echo "FAIL: ${dir##*/} constructed no mechanical constitutive law"
        return 1
    fi

    # The input must have come from the temperature case, and at the last time
    # too, not merely once at the start. OpenFOAM.com quotes the path it logs
    # and the other forks do not, so the quote is optional; the end of the
    # line is not, or time 4 would match 40
    if ! grep "${FRAMEWORK_T_READ}" "${dir}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | grep -Eq "hotCylinderTemperatureField/${COMPARISON_END_TIME}\"?$"
    then
        echo "FAIL: ${dir##*/} did not read T from" \
            "hotCylinderTemperatureField/${COMPARISON_END_TIME}"
        return 1
    fi

    local t
    t=$(solids4Foam::latestTime "${dir}")

    if ! awk "BEGIN {exit !((${t:-0} - ${COMPARISON_END_TIME})^2 <= 1e-20)}"
    then
        echo "FAIL: ${dir##*/} stopped at '${t}';" \
            "expected ${COMPARISON_END_TIME}"
        return 1
    fi

    local sigma
    sigma=$(extract_last "Max sigmaEq (von Mises stress)" "${dir}")

    if [[ -z "${sigma}" ]] \
      || ! awk "BEGIN {exit !(${sigma} >= ${SIGMA_NONZERO_MIN})}"
    then
        echo "FAIL: ${dir##*/} has no thermal stress (max sigmaEq '${sigma}')"
        return 1
    fi
}

# ------------------------------------------------------------
# Clean & run the arms
# ------------------------------------------------------------

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

SINGLE_DIR="${REGRESSION_ROOT}/single"

if [ "$CHECK_ONLY" = false ]; then
    for dir in "${CASE_DIR}" "${SINGLE_DIR}"; do
        copy_case "${dir}"
    done

    use_single_material "${SINGLE_DIR}"

    for dir in "${CASE_DIR}" "${SINGLE_DIR}"; do
        ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check_arm "${CASE_DIR}" || failures=$((failures + 1))
check_arm "${SINGLE_DIR}" || failures=$((failures + 1))

epsilon=$(extract_last "Max epsilonEq" "${CASE_DIR}" || true)
sigma=$(extract_last "Max sigmaEq (von Mises stress)" "${CASE_DIR}" || true)

if [[ -z "${epsilon}" || -z "${sigma}" ]]; then
    echo "FAIL: Could not extract one or more regression quantities"
    failures=$((failures + 1))
else
    if awk "BEGIN {exit !(${epsilon} >= ${EPS_MIN} && ${epsilon} <= ${EPS_MAX})}"
    then
        printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
    else
        printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${sigma} <= ${SIGMA_MAX})}"
    then
        printf "PASS: Max sigmaEq = %.6g\n" "${sigma}"
    else
        printf "FAIL: Max sigmaEq = %.6g\n" "${sigma}"
        failures=$((failures + 1))
    fi
fi

solids4Foam::checkFieldNorms "two materials, D" \
    "${CASE_DIR}/${COMPARISON_END_TIME}/D" \
    "${REF_D_MAX}" "${REF_D_MEAN}" "${FRAMEWORK_D_REL_TOL}" \
    || failures=$((failures + 1))

solids4Foam::checkFieldNorms "one material, D" \
    "${SINGLE_DIR}/${COMPARISON_END_TIME}/D" \
    "${REF_SINGLE_D_MAX}" "${REF_SINGLE_D_MEAN}" "${SINGLE_D_REL_TOL}" \
    || failures=$((failures + 1))

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
    for dir in "${CASE_DIR}" "${SINGLE_DIR}"; do
        ( cd "${dir}" && ./Allclean > /dev/null 2>&1 ) || true
    done
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
