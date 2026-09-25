#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# longWall regression test
# Uses top-surface stress and displacement histories.
#
# The implicit stiffness comes from the mechanicalConstitutiveLaw framework.
# impK is the coefficient of a Laplacian that is added implicitly and
# subtracted explicitly, so it sets how the solution is reached and not what it
# is. The answer must therefore be the one the removed legacy mechanicalModel
# reached, and that agreement is the check on the framework's finite-strain
# scalar tangent for MooneyRivlinElastic.
#
# The case also carries the framework's own checks, since its law is
# finite-strain only.
# ============================================================

UY_MIN=0.405
UY_MAX=0.407
SYY_MIN=9.99e7
SYY_MAX=1.001e8

# The final top-surface uy and sigma_yy of the removed legacy mechanicalModel,
# from the last commit that had it (mcl-stage8-coverage, c3a92b3d), the same on
# every fork.
#
# This case used to set solvePressureEqn, which made the legacy
# MooneyRivlinElastic solve a Laplacian equation for its hydrostatic stress.
# That smoothing stabilises the discretisation rather than describing the
# material, and with it the converged answers of the two models agreed to
# about 2e-6 relative, which is the evidence that it was only a stabilisation.
# The case does not need it, so it no longer sets it, and without it the legacy
# run matched the framework one in every written digit of D. The tolerance is
# left at the 1e-4 that accommodated the smoothing
LEGACY_UY=0.405906
LEGACY_SYY=1e+08
LEGACY_REL_TOL=1e-4

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
CONSTITUTIVE_LOGFILE="log.Test-mechanicalConstitutiveLaw"

APPROACHES=(
    main
)

echo "============================================================"
echo "longWall regression test"
echo "Top-surface uy in [${UY_MIN}, ${UY_MAX}] m"
echo "Top-surface sigma_yy in [${SYY_MIN}, ${SYY_MAX}] Pa"
echo "uy and sigma_yy match the legacy model to ${LEGACY_REL_TOL} relative"
echo "============================================================"
echo

prepare_case() {
    local approach="$1"
    local case_dir="${REGRESSION_ROOT}/${approach}"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done
}

find_history_file() {
    local case_dir="$1"
    local name="$2"

    find "${case_dir}/postProcessing" -name "${name}" -print 2>/dev/null \
        | tail -n 1
}

# The framework's own checks. This case's law is MooneyRivlinElastic, which
# implements no small-strain evaluation, so this is its runtime coverage
run_constitutive_test() {
    local case_dir="$1"

    # A skip where the application is not built, and a failure in CI, where
    # it always is
    solids4Foam::requireTestApp Test-mechanicalConstitutiveLaw \
        || return $(( $? - 1 ))

    if [[ ! -d "${case_dir}/constant/polyMesh" ]]; then
        echo "SKIP: mechanicalConstitutiveLaw checks (case has no mesh)"
        return 0
    fi

    if ( cd "${case_dir}" && Test-mechanicalConstitutiveLaw \
            > "${CONSTITUTIVE_LOGFILE}" 2>&1 )
    then
        local n_passed
        n_passed=$(grep -c 'PASS:' "${case_dir}/${CONSTITUTIVE_LOGFILE}" || true)

        if (( n_passed == 0 )); then
            echo "FAIL: mechanicalConstitutiveLaw checks reported no checks"
            return 1
        fi

        echo "PASS: mechanicalConstitutiveLaw checks (${n_passed} checks)"
        return 0
    fi

    echo "FAIL: mechanicalConstitutiveLaw checks"
    grep 'FAIL:' "${case_dir}/${CONSTITUTIVE_LOGFILE}" || true
    return 1
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

failures=0
constitutive_tested=false
declare -A RESULT_UY
declare -A RESULT_SYY

for approach in "${APPROACHES[@]}"; do
    CASE_DIR="${REGRESSION_ROOT}/${approach}"

    echo
    echo "------------------------------------------------------------"
    echo "Testing approach: ${approach}"
    echo "------------------------------------------------------------"

    if [ "$CHECK_ONLY" = false ]; then
        prepare_case "${approach}"
        ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
    else
        echo "Running in check-only mode: skipping Allclean and Allrun"
    fi

    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping ${approach} because it is unavailable in this environment"
        continue
    fi

    # impK must have come from the framework
    marker='Implicit stiffness from the mechanicalConstitutiveLaw framework'
    if grep -q "${marker}" "${CASE_DIR}/${SOLVER_LOGFILE}"; then
        echo "PASS: ${approach} took impK from the framework"
    else
        echo "FAIL: ${approach} did not take impK from the framework"
        failures=$((failures + 1))
    fi

    disp_file=$(find_history_file "${CASE_DIR}" 'solidDisplacementstop.dat')
    stress_file=$(find_history_file "${CASE_DIR}" 'solidStressestop.dat')

    if [[ -z "${disp_file}" || -z "${stress_file}" ]]; then
        echo "FAIL: ${approach} could not find one or more history files"
        failures=$((failures + 1))
        continue
    fi

    top_uy=$(awk 'END {print $9}' "${disp_file}")
    top_syy=$(awk 'END {print $5}' "${stress_file}")

    if [[ -z "${top_uy}" || -z "${top_syy}" ]]; then
        echo "FAIL: ${approach} could not extract top-surface values"
        failures=$((failures + 1))
        continue
    fi

    RESULT_UY["${approach}"]="${top_uy}"
    RESULT_SYY["${approach}"]="${top_syy}"

    if awk "BEGIN {exit !(${top_uy} >= ${UY_MIN} && ${top_uy} <= ${UY_MAX})}"; then
        printf "PASS: %s top-surface uy = %.6g\n" "${approach}" "${top_uy}"
    else
        printf "FAIL: %s top-surface uy = %.6g\n" "${approach}" "${top_uy}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${top_syy} >= ${SYY_MIN} && ${top_syy} <= ${SYY_MAX})}"; then
        printf "PASS: %s top-surface sigma_yy = %.6g\n" "${approach}" "${top_syy}"
    else
        printf "FAIL: %s top-surface sigma_yy = %.6g\n" "${approach}" "${top_syy}"
        failures=$((failures + 1))
    fi

    if [[ "${constitutive_tested}" == false ]]; then
        constitutive_tested=true

        if ! run_constitutive_test "${CASE_DIR}"; then
            failures=$((failures + 1))
        fi
    fi

    if [ "$CHECK_ONLY" = false ]; then
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    fi
done

# impK changes the path, not the answer
if [[ -n "${RESULT_UY[main]:-}" ]]; then
    for quantity in uy syy; do
        if [[ "${quantity}" == "uy" ]]; then
            a="${LEGACY_UY}"
            b="${RESULT_UY[main]}"
        else
            a="${LEGACY_SYY}"
            b="${RESULT_SYY[main]}"
        fi

        if awk "BEGIN {exit !(($a - $b)^2 <= (${LEGACY_REL_TOL}*$a)^2)}"; then
            printf "PASS: %s matches the legacy model (%.8g vs %.8g)\n" \
                "${quantity}" "$b" "$a"
        else
            printf "FAIL: %s differs from the legacy model (%.8g vs %.8g)\n" \
                "${quantity}" "$b" "$a"
            failures=$((failures + 1))
        fi
    done
elif [ "$CHECK_ONLY" = false ] \
    && ! solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"
then
    echo "FAIL: the main arm produced nothing to compare with the legacy model"
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
