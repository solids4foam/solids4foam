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

# GNU sed, for the in-place edits below
solids4Foam::requireGnuSed

# ============================================================
# suctionCaission regression test
#
# Pullout of a suction caisson: poroMechanicalLaw over
# linearElasticMohrCoulombPlastic, on poroLinearGeometry.
#
# The bands are wide, and deliberately so. This case does not
# converge tightly - it reaches the corrector limit on every
# time step, as it did with the removed legacy mechanicalModel
# - so the bands say the caisson yielded and suction developed,
# not that a particular number came back.
# ============================================================

# With the leastSquares gradient on every fork: a max sigmaEq of 1.05e6
# on OpenFOAM.com v2412, 9.0e5 on OpenFOAM.org 9 and 7.3e5 on foam-extend 4.1
EPS_MIN=0.85
EPS_MAX=1.25
SIGMA_MIN=6e5
SIGMA_MAX=2.2e6

# The comparison with the removed legacy mechanicalModel runs to t = 1 rather
# than ten. The full case takes about nine minutes, and the plastic return
# mapping and the effective stress the composite carries are both exercised
# within the first step
COMPARISON_END_TIME=1

# D and porePressure at t = 1, as the max and mean component magnitude of each
# field written to fourteen figures, per fork, with the leastSquares gradient
# on every fork: before, OpenFOAM ran pointCellsLeastSquares, which
# gave a max porePressure of 94942, and foam-extend leastSquares. The forks
# still differ, by 27% in the max porePressure, so the gradient scheme is not
# the cause of that difference; this case reaches the corrector limit
# on every step, as noted above. These are recorded numbers, though, and
# another compiler, CPU or MPI build moves an iterative solution by round-off
# at the solver tolerance: CI measures up to 3e-8 relative against values
# recorded on macOS. The tolerance, 1e-6 of the largest value, allows for
# that, and is ten times below the 1e-5 that a 0.001% change in a material
# constant makes
case "$(solids4Foam::foamFlavour)" in
    com)
        REF_D_MAX=0.039971853415323
        REF_D_MEAN=0.00326804246128935
        REF_P_MAX=91160.280037436
        REF_P_MEAN=19082.323226845
        ;;
    org)
        REF_D_MAX=0.039971831329265
        REF_D_MEAN=0.00326747514301296
        REF_P_MAX=91153.098748193
        REF_P_MEAN=19081.5856382194
        ;;
    foamextend)
        REF_D_MAX=0.039984257366288
        REF_D_MEAN=0.00342428784043281
        REF_P_MAX=66277.760264684
        REF_P_MEAN=17317.8068986057
        ;;
esac
REF_REL_TOL=1e-6

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "suctionCaission regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "Plus the comparison with the reference, to t=${COMPARISON_END_TIME}"
echo "============================================================"
echo

prepare_case() {
    local dir="$1"

    rm -rf "${dir}"
    mkdir -p "${dir}"

    local item base_item
    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${dir}/"
    done
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
    prepare_case "${CASE_DIR}"
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped here"
    exit 0
fi

# Run the case to the comparison time and hold it to the legacy answer there
run_reference_comparison() {
    local dir="${REGRESSION_ROOT}/comparison"

    prepare_case "${dir}"

    "${SOLIDS4FOAM_SED}" -i "s|^endTime .*|endTime         ${COMPARISON_END_TIME};|" \
        "${dir}/system/controlDict"

    # Enough digits that the comparison is about the solution and not about
    # the last figure written
    if grep -q "^writePrecision" "${dir}/system/controlDict"; then
        "${SOLIDS4FOAM_SED}" -i 's|^writePrecision.*|writePrecision  14;|' \
            "${dir}/system/controlDict"
    else
        echo "writePrecision  14;" >> "${dir}/system/controlDict"
    fi

    ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
        echo "FAIL: the comparison could not run ${dir}"
        return 1
    }

    if ! grep -q "Selecting mechanical constitutive law" \
        "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the comparison constructed no mechanical constitutive law"
        return 1
    fi

    local t
    t=$(solids4Foam::latestTime "${dir}")

    if [[ -z "${t}" ]] \
        || ! awk "BEGIN {exit !((${t} - ${COMPARISON_END_TIME})^2 <= 1e-20)}"
    then
        echo "FAIL: the comparison stopped at '${t}', not at ${COMPARISON_END_TIME}"
        return 1
    fi

    local failed=0

    solids4Foam::checkFieldNorms "D at t = ${t}" "${dir}/${t}/D" \
        "${REF_D_MAX}" "${REF_D_MEAN}" "${REF_REL_TOL}" \
        || failed=1

    solids4Foam::checkFieldNorms "porePressure at t = ${t}" \
        "${dir}/${t}/porePressure" \
        "${REF_P_MAX}" "${REF_P_MEAN}" "${REF_REL_TOL}" \
        || failed=1

    return "${failed}"
}

epsilon=$(grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null \
    | tail -n 1 | awk '{print $NF}' || true)
sigma=$(grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
    2>/dev/null | tail -n 1 | awk '{print $NF}' || true)

if [[ -z "${epsilon}" || -z "${sigma}" ]]; then
    echo "FAIL: could not extract the regression quantities"
    exit 1
fi

failures=0

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

if [ "$CHECK_ONLY" = false ] && ! run_reference_comparison; then
    failures=$((failures + 1))
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
