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
# suctionCaission regression test
#
# Pullout of a suction caisson: poroMechanicalLaw over
# linearElasticMohrCoulombPlastic, on poroLinearGeometry.
#
# The bands are wide, and deliberately so. This case does not
# converge tightly - it reaches the corrector limit on every
# time step, in the legacy solver as much as the framework one
# - so the bands say the caisson yielded and suction developed,
# not that a particular number came back.
# ============================================================

EPS_MIN=0.85
EPS_MAX=1.25
SIGMA_MIN=1.5e6
SIGMA_MAX=2.2e6

# The framework comparison runs two steps rather than ten. The full case takes
# about nine minutes an arm, and the plastic return mapping and the effective
# stress the composite carries are both exercised within the first step; ten
# steps of agreement would cost twenty minutes to say the same thing
COMPARISON_END_TIME=1

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "suctionCaission regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "Plus the legacy-versus-framework comparison, to t=${COMPARISON_END_TIME}"
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

# Run the case both ways and require the two to agree exactly. Both arms are
# shortened, and both are shortened the same way, so the comparison is still
# between two runs that differ in one dictionary entry and nothing else
run_framework_comparison() {
    local legacy_dir="${REGRESSION_ROOT}/comparisonLegacy"
    local framework_dir="${REGRESSION_ROOT}/comparisonFramework"
    local dir

    for dir in "${legacy_dir}" "${framework_dir}"; do
        prepare_case "${dir}"

        sed -i "s|^endTime .*|endTime         ${COMPARISON_END_TIME};|" \
            "${dir}/system/controlDict"

        # Enough digits that the comparison is about the solution and not
        # about the last figure written
        if grep -q "^writePrecision" "${dir}/system/controlDict"; then
            sed -i 's|^writePrecision.*|writePrecision  14;|' \
                "${dir}/system/controlDict"
        else
            echo "writePrecision  14;" >> "${dir}/system/controlDict"
        fi
    done

    sed -i \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${framework_dir}/constant/solidProperties"

    for dir in "${legacy_dir}" "${framework_dir}"; do
        ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
            echo "FAIL: the comparison could not run ${dir}"
            return 1
        }
    done

    # Each arm must have taken the path it was set up for
    if ! grep -q "Selecting mechanical constitutive law" \
        "${framework_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the framework arm did not use the framework"
        return 1
    fi

    if grep -q "Selecting mechanical constitutive law" \
        "${legacy_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the legacy arm used the framework"
        return 1
    fi

    local tL tF
    tL=$(foamListTimes -case "${legacy_dir}" -latestTime 2>/dev/null | tail -n 1)
    tF=$(foamListTimes -case "${framework_dir}" -latestTime 2>/dev/null \
        | tail -n 1)

    if [[ -z "${tL}" || "${tL}" != "${tF}" ]]; then
        echo "FAIL: the arms reached different times ('${tL}' vs '${tF}')"
        return 1
    fi

    local ok=0
    local f
    for f in D p; do
        if [[ ! -f "${legacy_dir}/${tL}/${f}" ]]; then
            continue
        fi

        if ! diff -q "${legacy_dir}/${tL}/${f}" "${framework_dir}/${tF}/${f}" \
            > /dev/null
        then
            echo "FAIL: framework and legacy differ in ${f}"
            return 1
        fi

        ok=$((ok + 1))
    done

    if (( ok == 0 )); then
        echo "FAIL: the comparison produced no fields"
        return 1
    fi

    echo "PASS: framework and legacy agree exactly (${ok} fields)"
    return 0
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

if [ "$CHECK_ONLY" = false ] && ! run_framework_comparison; then
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
