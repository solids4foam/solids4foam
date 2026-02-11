#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

# -----------------------------------------------------------------------------
# Regression test for rigid rotation of a hyperelastic sphere
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

# Solution approaches to test
APPROACHES=(
    totalLagrangian
    updatedLagrangian
    totalLagrangianPetscSnes
    updatedLagrangianPetscSnes
)

failures=0

echo "============================================================"
echo "Rigid rotation regression test"
echo "Stress threshold: sigmaEq < ${SIGMA_TOL}"
echo "============================================================"

for approach in "${APPROACHES[@]}"; do
    echo
    echo "------------------------------------------------------------"
    echo "Testing approach: ${approach}"
    echo "------------------------------------------------------------"

    # Clean previous run
    ( ./Allclean ) >/dev/null 2>&1 || true

    # Run case
    ( ./Allrun "${approach}" ) > "${ALLRUN_LOGFILE}" 2>&1

    # Extract final Max sigmaEq
    sigma=$(grep "Max sigmaEq (von Mises stress)" "${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1 || true)

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
( ./Allclean ) >/dev/null 2>&1 || true

echo
echo "============================================================"

if (( failures > 0 )); then
    echo "Regression test FAILED (${failures} failing case(s))"
    exit 1
else
    echo "Regression test PASSED (all approaches)"
    exit 0
fi
