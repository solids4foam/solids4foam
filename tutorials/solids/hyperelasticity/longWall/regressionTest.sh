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
# Run twice: once with the implicit stiffness from the legacy mechanicalModel
# and once with it from the mechanicalConstitutiveLaw framework. impK is the
# coefficient of a Laplacian that is added implicitly and subtracted
# explicitly, so it sets how the solution is reached and not what it is. The
# two runs must therefore agree, and that agreement is the check on the
# framework's finite-strain scalar tangent for MooneyRivlinElastic.
#
# The case also carries the framework's own checks, since its law is
# finite-strain only.
# ============================================================

UY_MIN=0.405
UY_MAX=0.407
SYY_MIN=9.99e7
SYY_MAX=1.001e8

# The two runs are NOT expected to be bit-identical here, unlike
# rotatingCylinder. This case sets solvePressureEqn, so the legacy
# MooneyRivlinElastic solves a Laplacian equation for its hydrostatic stress,
# which the framework law deliberately omits: that smoothing stabilises the
# discretisation rather than describing the material, and belongs to the solid
# model. The converged answers agree to about 2e-6 relative, which is the
# evidence that it is indeed a stabilisation, so the tolerance is set to
# accommodate that rather than to hide it
CROSS_TOL=1e-4

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
CONSTITUTIVE_LOGFILE="log.Test-mechanicalConstitutiveLaw"

APPROACHES=(
    legacy
    framework
)

echo "============================================================"
echo "longWall regression test"
echo "Top-surface uy in [${UY_MIN}, ${UY_MAX}] m"
echo "Top-surface sigma_yy in [${SYY_MIN}, ${SYY_MAX}] Pa"
echo "Legacy and framework impK agree to ${CROSS_TOL} relative"
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

    if [[ "${approach}" == "framework" ]]; then
        # The switch is read from solidModelDict(), which is the <type>Coeffs
        # sub-dictionary, so it must go inside those braces. Appended at the
        # top level it is silently ignored and this approach would quietly be
        # a second legacy run
        sed -i.bak \
            's/^    nCorrectors\(.*\)$/    useMechanicalConstitutiveLawManager yes;\n    nCorrectors\1/' \
            "${case_dir}/constant/solidProperties"
        rm -f "${case_dir}/constant/solidProperties.bak"

        if ! grep -q 'useMechanicalConstitutiveLawManager' \
            "${case_dir}/constant/solidProperties"
        then
            echo "FAIL: could not enable the framework in solidProperties"
            exit 1
        fi
    fi
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

    if ! command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
        echo "SKIP: Test-mechanicalConstitutiveLaw not found in PATH"
        return 0
    fi

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
            echo "SKIP: mechanicalConstitutiveLaw checks (no checks reported)"
            return 0
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

    # Assert the run really took the path this approach names
    marker='Implicit stiffness from the mechanicalConstitutiveLaw framework'
    if grep -q "${marker}" "${CASE_DIR}/${SOLVER_LOGFILE}"; then
        used_framework=true
    else
        used_framework=false
    fi

    if [[ "${approach}" == "framework" && "${used_framework}" == false ]]; then
        echo "FAIL: framework approach did not use the framework impK"
        failures=$((failures + 1))
    elif [[ "${approach}" == "legacy" && "${used_framework}" == true ]]; then
        echo "FAIL: legacy approach unexpectedly used the framework impK"
        failures=$((failures + 1))
    else
        echo "PASS: ${approach} took the expected impK path"
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

# The point of running both: impK changes the path, not the answer
if [[ -n "${RESULT_UY[legacy]:-}" && -n "${RESULT_UY[framework]:-}" ]]; then
    for quantity in uy syy; do
        if [[ "${quantity}" == "uy" ]]; then
            a="${RESULT_UY[legacy]}"
            b="${RESULT_UY[framework]}"
        else
            a="${RESULT_SYY[legacy]}"
            b="${RESULT_SYY[framework]}"
        fi

        if awk "BEGIN {exit !(($a - $b)^2 <= ($CROSS_TOL*$a)^2)}"; then
            printf "PASS: legacy and framework %s agree (%.8g vs %.8g)\n" \
                "${quantity}" "$a" "$b"
        else
            printf "FAIL: legacy and framework %s differ (%.8g vs %.8g)\n" \
                "${quantity}" "$a" "$b"
            failures=$((failures + 1))
        fi
    done
else
    echo "SKIP: cross-check needs both approaches to have run"
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
