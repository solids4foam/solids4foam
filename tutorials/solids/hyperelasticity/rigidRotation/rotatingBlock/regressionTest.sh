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

# The high-order arms take their stress at the face quadrature points from the
# mechanicalConstitutiveLaw framework. Their final D, as written, matched the
# removed legacy mechanicalModel's exactly, and still must. It is recorded as
# the max and mean component magnitude of its internal values, from the last
# commit that had the legacy model (mcl-stage8-coverage, c3a92b3d), on
# OpenFOAM.com v2512 and OpenFOAM.org 9, where the high-order arms run, and the
# same on both. The tolerance, 1e-5 of the largest value, is a few units in the
# last of the six figures written.
#
# At those six figures, D is the rigid rotation itself: every arm here writes
# the same field. So this says the high-order arms still rotate the block
# rigidly to that precision, as they did on the legacy model; the stress bound
# above is the check on how rigidly
declare -A REF_D_MAX=(
    [highOrder]=2.05053
    [highOrderUpdatedLagrangian]=2.05053
)
declare -A REF_D_MEAN=(
    [highOrder]=0.4825039126
    [highOrderUpdatedLagrangian]=0.4825039126
)
REF_D_REL_TOL=1e-5

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

    if ! grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${approach} constructed no mechanical constitutive law"
        failures=$((failures + 1))
    fi

    if [[ -n "${REF_D_MAX[${approach}]:-}" ]]; then
        latest_time=$(solids4Foam::latestTime "${CASE_DIR}")

        if ! solids4Foam::checkFieldNorms "${approach} D" \
            "${CASE_DIR}/${latest_time}/D" "${REF_D_MAX[${approach}]}" \
            "${REF_D_MEAN[${approach}]}" "${REF_D_REL_TOL}"
        then
            failures=$((failures + 1))
        fi
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
