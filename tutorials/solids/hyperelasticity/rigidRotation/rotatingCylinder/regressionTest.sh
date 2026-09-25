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
# The implicit stiffness comes from the mechanicalConstitutiveLaw framework.
# impK is the coefficient of the Laplacian term, so it affects how the solution
# is reached and not what it is, and the answer must be the one the removed
# legacy mechanicalModel reached.
#
# The case also carries the framework's own checks, because its law is
# StVenantKirchhoffElastic, which is finite-strain only.
# -----------------------------------------------------------------------------

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
CONSTITUTIVE_LOGFILE="log.Test-mechanicalConstitutiveLaw"

# Stress threshold (deliberately loose)
SIGMA_TOL=1e4

# The final sigmaEq of the removed legacy mechanicalModel, from the last commit
# that had it (mcl-stage8-coverage, c3a92b3d), the same on every fork. It is
# near zero, so it is compared on the absolute scale of the threshold - 1% of
# it - rather than relative to a value that is itself almost zero, as the two
# models were
LEGACY_SIGMA=3700.73
LEGACY_SIGMA_TOL=$(awk "BEGIN {print 0.01*${SIGMA_TOL}}")

APPROACHES=(
    main
)

failures=0
declare -A RESULT_SIGMA

echo "============================================================"
echo "Rigid rotation cylinder regression test"
echo "Stress threshold: sigmaEq < ${SIGMA_TOL}"
echo "sigmaEq within ${LEGACY_SIGMA_TOL} of the legacy model"
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
}

# Exercise the framework's own checks on this case. Its law is
# StVenantKirchhoffElastic, which implements no small-strain evaluation, so
# this is the runtime coverage of that law's finite-strain paths
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

    marker='Implicit stiffness from the mechanicalConstitutiveLaw framework'
    if grep -q "${marker}" "${CASE_DIR}/${SOLVER_LOGFILE}"; then
        echo "PASS: ${approach} took impK from the framework"
    else
        echo "FAIL: ${approach} did not take impK from the framework"
        failures=$((failures + 1))
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

# impK changes the path, not the answer
if [[ -n "${RESULT_SIGMA[main]:-}" ]]
then
    a="${LEGACY_SIGMA}"
    b="${RESULT_SIGMA[main]}"

    if awk "BEGIN {exit !(($a - $b)^2 < (${LEGACY_SIGMA_TOL})^2)}"; then
        printf "PASS: sigmaEq matches the legacy model (%.6g vs %.6g)\n" \
            "$b" "$a"
    else
        printf "FAIL: sigmaEq differs from the legacy model (%.6g vs %.6g)\n" \
            "$b" "$a"
        failures=$((failures + 1))
    fi
elif ! solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"
then
    echo "FAIL: the main arm produced no sigmaEq to compare with the legacy model"
    failures=$((failures + 1))
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
