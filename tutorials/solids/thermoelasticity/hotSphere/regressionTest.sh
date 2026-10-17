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
# hotSphere regression test
# Uses the tutorial's reported temperature and stress extrema.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

T_MIN=339.5
T_MAX=340.5
EPS_MIN=1.5e-4
EPS_MAX=2.5e-4
SIGMA_MIN=3.8e7
SIGMA_MAX=4.6e7

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "hotSphere regression test"
echo "Max T       in [${T_MIN}, ${T_MAX}]"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq in [${SIGMA_MIN}, ${SIGMA_MAX}]"
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

shorten_case() {
    local controlDict="${CASE_DIR}/system/controlDict"
    "${SOLIDS4FOAM_SED}" -i.bak 's/^endTime[[:space:]]\+5;/endTime         1;/' "${controlDict}"
    rm -f "${controlDict}.bak"
}

# ------------------------------------------------------------
# Clean & run case
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

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    shorten_case
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_temperature() {
    grep "Max T" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk '{print $NF}'
}

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk '{print $NF}'
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk '{print $NF}'
}

# Run the case to its full end time, and hold it to the answer of the removed
# legacy mechanicalModel.
#
# This case runs several time steps, so unlike slabCooling it exercises the
# framework rolling its constitutive state over between them.
#
# The comparison is to a tolerance rather than exact. The framework and the
# legacy model solved the same problem but reached it by slightly different
# iteration paths, because the implicit stiffness that steers the iteration is
# built differently: the framework interpolates its cell tangent to the faces
# where the legacy model formed a face value directly. The converged answers
# therefore agreed only to the solution tolerance, and the difference shrank
# with it - at solutionTolerance 1e-6 it was around 1e-7 of the displacement,
# and tightening to 1e-10 took it to 1e-8. The threshold here is well above
# that and far below anything physical
FRAMEWORK_D_REL_TOL=1e-6
COMPARISON_END_TIME=5

# The final D, as the max and mean component magnitude of the field written
# to fourteen figures, from the removed legacy mechanicalModel on the last
# commit that had it (mcl-stage8-coverage, c3a92b3d). OpenFOAM.com v2512's;
# OpenFOAM.org 9 agrees to 2e-9. foam-extend is not compared, as below
REF_D_MAX=0.00015702498751748
REF_D_MEAN=7.86486280639836e-05

run_reference_comparison() {
    # Not on foam-extend, where the framework and the legacy model differed by
    # 0.6 % in D, and the framework's answer is the one kept. They matched to
    # 1e-13 for two correctors and part of the third, and differed only in the
    # stress on the three symmetryPlane patches. The framework corrects the
    # stress's boundary conditions after evaluating it; the legacy law did not.
    # On foam-extend's symmetryPlane that correction replaces the law's
    # boundary value with the symmetry transform of the adjacent cell's
    # stress, which has zero shear traction on the plane, as symmetry
    # requires. Correcting the legacy stress too reproduced the framework's
    # answer to 2e-8.
    #
    # The transform was checked because older OpenFOAM versions were thought to
    # get it wrong for tensors: 0.5*(x + transform(I - 2nn, x)) matches a
    # hand-written reference to 4e-16 for symmTensor and tensor, with zero shear
    # traction and all other components kept, on v1912-v2606, OpenFOAM.org 8-13
    # and foam-extend 4.1 and 5.0. So the legacy answer there is the one known
    # to be wrong, and there is nothing to hold the framework to
    if [[ "${WM_PROJECT:-}" == "foam" ]]; then
        echo "SKIP: reference comparison (known foam-extend symmetryPlane difference)"
        return 0
    fi

    local dir="${REGRESSION_ROOT}/comparison"

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

    # Write enough digits for the comparison to be about the solution rather
    # than about the file format. At the default six significant figures the
    # two differ by around 3e-6 simply because that is the last digit written,
    # which would tell us nothing
    if grep -q "^writePrecision" "${dir}/system/controlDict"; then
        "${SOLIDS4FOAM_SED}" -i 's|^writePrecision.*|writePrecision  14;|' \
            "${dir}/system/controlDict"
    else
        echo "writePrecision  14;" >> "${dir}/system/controlDict"
    fi

    ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
        echo "FAIL: the reference comparison could not run ${dir}"
        return 1
    }

    if ! grep -q "Selecting mechanical constitutive law" \
        "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the comparison constructed no mechanical constitutive law"
        return 1
    fi

    if ! grep -q "^End" "${dir}/${SOLVER_LOGFILE}" \
      || grep -qE "Nonlinear solve did not converge|SNES convergence error|FOAM FATAL" \
          "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${dir} did not complete and converge"
        return 1
    fi

    local t
    t=$(solids4Foam::latestTime "${dir}")

    if [[ -z "${t}" ]] \
        || ! awk "BEGIN {exit !((${t} - ${COMPARISON_END_TIME})^2 <= 1e-20)}"
    then
        echo "FAIL: comparison stopped at '${t}'; expected ${COMPARISON_END_TIME}"
        return 1
    fi

    solids4Foam::checkFieldNorms "D at t = ${t}" "${dir}/${t}/D" \
        "${REF_D_MAX}" "${REF_D_MEAN}" "${FRAMEWORK_D_REL_TOL}"
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

temperature=$(extract_max_temperature)
epsilon=$(extract_max_epsilon)
sigma=$(extract_max_sigma)

if [[ -z "${temperature}" || -z "${epsilon}" || -z "${sigma}" ]]; then
    echo "FAIL: Could not extract one or more regression quantities"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${temperature} >= ${T_MIN} && ${temperature} <= ${T_MAX})}"; then
    printf "PASS: Max T = %.6g\n" "${temperature}"
else
    printf "FAIL: Max T = %.6g\n" "${temperature}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${epsilon} >= ${EPS_MIN} && ${epsilon} <= ${EPS_MAX})}"; then
    printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
else
    printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${SIGMA_MAX} >= ${sigma})}"; then
    printf "PASS: Max sigmaEq = %.6g\n" "${sigma}"
else
    printf "FAIL: Max sigmaEq = %.6g\n" "${sigma}"
    failures=$((failures + 1))
fi

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
    if ! run_reference_comparison; then
        failures=$((failures + 1))
    fi

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
