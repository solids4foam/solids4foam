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
# removed legacy mechanicalModel's exactly, and still must. The legacy fields
# are recorded as digests of their written internal values, from the last
# commit that had the legacy model (mcl-stage8-coverage, c3a92b3d), on
# OpenFOAM.com v2512 and OpenFOAM.org 9, where the high-order arms run, and the
# same on both.
#
# At the six figures written, D is the rigid rotation itself: every arm here
# writes the same field. So this says the high-order arms still rotate the
# block rigidly to that precision, as they did on the legacy model; the
# stress bound above is the check on how rigidly
declare -A LEGACY_D_DIGEST=()
case "$(solids4Foam::foamFlavour)" in
    com|org)
        LEGACY_D_DIGEST[highOrder]=59bbcbdcd22ae9b18f041bbbc1d5916cefeb9adb
        LEGACY_D_DIGEST[highOrderUpdatedLagrangian]=59bbcbdcd22ae9b18f041bbbc1d5916cefeb9adb
        ;;
esac

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

# A digest of a field's internal values as written, independent of the file
# header, which names the OpenFOAM version that wrote it
internal_field_digest() {
    python3 - "$1" << 'PYEOF'
import hashlib
import re
import sys

text = open(sys.argv[1]).read()
field = re.search(r"\binternalField\s+(.*?);\s*\n\s*boundaryField", text, re.DOTALL)
if not field:
    sys.exit(f"cannot find the internalField in {sys.argv[1]}")
print(hashlib.sha1(" ".join(field.group(1).split()).encode()).hexdigest())
PYEOF
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

    if [[ -n "${LEGACY_D_DIGEST[${approach}]:-}" ]]; then
        latest_time=$(solids4Foam::latestTime "${CASE_DIR}")
        digest=$(internal_field_digest "${CASE_DIR}/${latest_time}/D" || true)

        if [[ "${digest}" == "${LEGACY_D_DIGEST[${approach}]}" ]]; then
            echo "PASS: ${approach} D matches the legacy model's exactly"
        else
            echo "FAIL: ${approach} D differs from the legacy model's" \
                "(digest ${digest:-none}, legacy ${LEGACY_D_DIGEST[${approach}]})"
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
