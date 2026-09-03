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
# stripFooting regression test
#
# Bearing capacity of a strip footing, poroMechanicalLaw over
# linearElasticMohrCoulombPlastic. The bands below are wide
# enough to be about the physics rather than about round-off:
# the case must reach failure, so the strain has to be well
# into the plastic range, and the excess pore pressure has to
# be of the right order.
# ============================================================

EPS_MIN=0.012
EPS_MAX=0.021
P_MIN=5.0e4
P_MAX=9.0e4

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "stripFooting regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max |p|       in [${P_MIN}, ${P_MAX}]"
echo "Plus the legacy-versus-framework comparison"
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

    # Enough digits that a comparison is about the solution rather than about
    # the last figure written
    if grep -q "^writePrecision" "${dir}/system/controlDict"; then
        sed -i 's|^writePrecision.*|writePrecision  14;|' \
            "${dir}/system/controlDict"
    else
        echo "writePrecision  14;" >> "${dir}/system/controlDict"
    fi
}

max_abs_field() {
    # Largest magnitude of any component of a field's internal values
    awk '
        function abs(x) { return x < 0 ? -x : x }
        /^internalField/ { inField = 1 }
        inField && /^\(/ { inList = 1; next }
        inList && /^\)/ { inList = 0 }
        inList {
            gsub(/[()]/, "")
            for (i = 1; i <= NF; i++)
            {
                v = abs($i); if (v > m) m = v
            }
        }
        END { printf "%.10g\n", m }
    ' "$1"
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

# Run the case again with the stress taken from the mechanicalConstitutiveLaw
# framework rather than the legacy mechanicalModel, and require the two to
# agree exactly.
#
# This is the case that exercises the Mohr-Coulomb return mapping and the
# history it carries, over thirty-eight steps, underneath the poro composite
run_framework_comparison() {
    local dir="${REGRESSION_ROOT}/framework"

    prepare_case "${dir}"

    # The two arms differ in this one entry and nothing else
    sed -i \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${dir}/constant/solidProperties"

    ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
        echo "FAIL: the framework arm did not run"
        return 1
    }

    # Each arm must have taken the path it was set up for
    if ! grep -q "Selecting mechanical constitutive law" \
        "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the framework arm did not use the framework"
        return 1
    fi

    if grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the legacy arm used the framework"
        return 1
    fi

    local tL tF
    tL=$(foamListTimes -case "${CASE_DIR}" -latestTime 2>/dev/null | tail -n 1)
    tF=$(foamListTimes -case "${dir}" -latestTime 2>/dev/null | tail -n 1)

    if [[ -z "${tL}" || "${tL}" != "${tF}" ]]; then
        echo "FAIL: the arms reached different times ('${tL}' vs '${tF}')"
        return 1
    fi

    if [[ ! -f "${CASE_DIR}/${tL}/D" || ! -f "${dir}/${tF}/D" ]]; then
        echo "FAIL: the comparison produced no D field"
        return 1
    fi

    if diff -q "${CASE_DIR}/${tL}/D" "${dir}/${tF}/D" > /dev/null; then
        echo "PASS: framework and legacy agree exactly"
        return 0
    fi

    echo "FAIL: framework and legacy differ"
    return 1
}

epsilon=$(grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null \
    | tail -n 1 | awk '{print $NF}' || true)

latest=$(foamListTimes -case "${CASE_DIR}" -latestTime 2>/dev/null | tail -n 1)

if [[ -z "${epsilon}" || -z "${latest}" ]]; then
    echo "FAIL: could not extract the regression quantities"
    exit 1
fi

pressure=$(max_abs_field "${CASE_DIR}/${latest}/p")

failures=0

if awk "BEGIN {exit !(${epsilon} >= ${EPS_MIN} && ${epsilon} <= ${EPS_MAX})}"
then
    printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
else
    printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${pressure} >= ${P_MIN} && ${pressure} <= ${P_MAX})}"
then
    printf "PASS: Max |p| = %.6g\n" "${pressure}"
else
    printf "FAIL: Max |p| = %.6g\n" "${pressure}"
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
