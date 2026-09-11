#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# -----------------------------------------------------------------------------
# Regression test for rigid rotation of a hyperelastic cylinder
#
# Physics invariant:
#   Pure rigid-body rotation should produce (near) zero stress.
#
# We check that the final reported Max sigmaEq remains below a loose threshold.
#
# The case is run twice: once with the implicit stiffness from the legacy
# mechanicalModel and once with it from the mechanicalConstitutiveLaw
# framework. impK is the coefficient of the Laplacian term, so it affects how
# the solution is reached and not what it is, and the two runs must agree.
#
# The case also carries the framework's own checks, because its law is
# StVenantKirchhoffElastic, which is finite-strain only.
# -----------------------------------------------------------------------------

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
CONSTITUTIVE_LOGFILE="log.Test-mechanicalConstitutiveLaw"

# Stress threshold (deliberately loose)
SIGMA_TOL=1e4

APPROACHES=(
    legacy
    framework
)

failures=0
declare -A RESULT_SIGMA

echo "============================================================"
echo "Rigid rotation cylinder regression test"
echo "Stress threshold: sigmaEq < ${SIGMA_TOL}"
echo "Legacy and framework impK both run and compared"
echo "============================================================"

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
            's/^    \/\/nCorrectors\(.*\)$/    useMechanicalConstitutiveLawManager yes;\n    \/\/nCorrectors\1/' \
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

# Exercise the framework's own checks on this case. Its law is
# StVenantKirchhoffElastic, which implements no small-strain evaluation, so
# this is the runtime coverage of that law's finite-strain paths
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

constitutive_tested=false

for approach in "${APPROACHES[@]}"; do
    CASE_DIR="${REGRESSION_ROOT}/${approach}"

    echo
    echo "------------------------------------------------------------"
    echo "Testing approach: ${approach}"
    echo "------------------------------------------------------------"

    prepare_case "${approach}"

    ( cd "${CASE_DIR}" && ./Allclean ) >/dev/null 2>&1 || true
    ( cd "${CASE_DIR}" && ./Allrun ) > "${CASE_DIR}/${ALLRUN_LOGFILE}" 2>&1

    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping ${approach} because it is unavailable in this environment"
        continue
    fi

    # Assert the run really took the path this approach names. Without this a
    # misplaced switch makes "framework" a second legacy run that passes
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

    sigma=$(grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1 || true)

    if [[ -z "${sigma}" ]]; then
        echo "FAIL: ${approach} could not extract sigmaEq from log"
        failures=$((failures + 1))
        continue
    fi

    RESULT_SIGMA["${approach}"]="${sigma}"

    if awk "BEGIN {exit !(${sigma} < ${SIGMA_TOL})}"; then
        printf "PASS: %s final sigmaEq = %.6g\n" "${approach}" "${sigma}"
    else
        printf "FAIL: %s final sigmaEq = %.6g exceeds threshold %.6g\n" \
            "${approach}" "${sigma}" "${SIGMA_TOL}"
        failures=$((failures + 1))
    fi

    # Before the Allclean below, which removes the mesh
    if [[ "${constitutive_tested}" == false ]]; then
        constitutive_tested=true

        if ! run_constitutive_test "${CASE_DIR}"; then
            failures=$((failures + 1))
        fi
    fi

    ( cd "${CASE_DIR}" && ./Allclean ) >/dev/null 2>&1 || true
done

# The point of running both: impK changes the path, not the answer. Both are
# near zero here, so compare on the absolute scale of the threshold rather
# than relative to a value that is itself almost zero
if [[ -n "${RESULT_SIGMA[legacy]:-}" && -n "${RESULT_SIGMA[framework]:-}" ]]
then
    a="${RESULT_SIGMA[legacy]}"
    b="${RESULT_SIGMA[framework]}"

    if awk "BEGIN {exit !(($a - $b)^2 < (0.01*$SIGMA_TOL)^2)}"; then
        printf "PASS: legacy and framework sigmaEq agree (%.6g vs %.6g)\n" \
            "$a" "$b"
    else
        printf "FAIL: legacy and framework sigmaEq differ (%.6g vs %.6g)\n" \
            "$a" "$b"
        failures=$((failures + 1))
    fi
else
    echo "SKIP: cross-check needs both approaches to have run"
fi

echo
echo "============================================================"

if (( failures > 0 )); then
    echo "Regression test FAILED"
    exit 1
else
    echo "Regression test PASSED"
    exit 0
fi
