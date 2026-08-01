#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# -----------------------------------------------------------------------------
# Regression test for rigid rotation of a hyperelastic block
#
# Physics invariant:
#   Pure rigid-body rotation should produce (near) zero stress.
#
# We check that the final reported Max sigmaEq remains below a loose threshold.
# -----------------------------------------------------------------------------

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

# Stress threshold (deliberately loose)
SIGMA_TOL=50.0

# Solution approach to test
APPROACHES=(
    totalLagrangian
    totalLagrangianPetscSnes
    updatedLagrangianPetscSnes
    highOrder
    highOrderUpdatedLagrangian
)

failures=0

echo "============================================================"
echo "Rigid rotation block regression test"
echo "Stress threshold: sigmaEq < ${SIGMA_TOL}"
echo "============================================================"

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

prepare_case

for approach in "${APPROACHES[@]}"; do
    echo
    echo "------------------------------------------------------------"
    echo "Testing approach: ${approach}"
    echo "------------------------------------------------------------"

    # Clean previous run
    # The solver log is removed explicitly so that a failed run cannot be
    # checked against the log of the previous approach
    ( cd "${CASE_DIR}" && ./Allclean ) >/dev/null 2>&1 || true
    rm -f "${CASE_DIR}/${SOLVER_LOGFILE}"

    # Run case
    ( cd "${CASE_DIR}" && ./Allrun "${approach}" ) > "${CASE_DIR}/${ALLRUN_LOGFILE}" 2>&1

    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        continue
    fi

    # Extract final Max sigmaEq
    sigma=$(grep "Max sigmaEq" "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null \
        | tail -n 1 \
        | awk '{print $NF}' || true)

    if [[ -z "${sigma}" ]]; then
        echo "FAIL: Could not extract sigmaEq from log"
        failures=$((failures + 1))
        continue
    fi

    # Compare using awk for floating-point safety
    if awk "BEGIN {exit !(${sigma} < ${SIGMA_TOL})}"; then
        printf "PASS: final sigmaEq = %.6g\n" "${sigma}"
    else
        printf "FAIL: final sigmaEq = %.6g exceeds threshold %.6g\n" \
            "${sigma}" "${SIGMA_TOL}"
        failures=$((failures + 1))
    fi
done

# Clean the case
( cd "${CASE_DIR}" && ./Allclean ) >/dev/null 2>&1 || true

echo
echo "============================================================"

if (( failures > 0 )); then
    echo "Regression test FAILED (${failures} failing case(s))"
    exit 1
else
    echo "Regression test PASSED (all approaches)"
    exit 0
fi
