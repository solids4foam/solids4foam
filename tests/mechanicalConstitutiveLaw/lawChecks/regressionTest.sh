#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# mechanicalConstitutiveLaw law checks
# Runs Test-mechanicalConstitutiveLawChecks on a one-cell mesh, which checks
# every law listed in constant/lawChecks against the properties a law of its
# kind must have: a stress-free reference state, a symmetric reference
# tangent, a scalar tangent equal to its normal stiffness where the law is
# isotropic, objectivity at finite strain and linearity where it is linear
# ============================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LOGFILE="log.Test-mechanicalConstitutiveLawChecks"

# Source solids4Foam scripts. From PATH, as a tutorial's Allrun does, since
# this case is run both in place and from a copy at a different depth
source solids4FoamScripts.sh

echo "============================================================"
echo "mechanicalConstitutiveLaw law checks"
echo "============================================================"

cd "${SCRIPT_DIR}"

# A skip where the application is not built, and a failure in CI, where it
# always is
solids4Foam::requireTestApp Test-mechanicalConstitutiveLawChecks \
    || exit $(( $? - 1 ))

rm -rf 0 constant/polyMesh "${LOGFILE}" log.blockMesh
mkdir -p 0

solids4Foam::convertCaseFormat .
blockMesh > log.blockMesh 2>&1

status=0
Test-mechanicalConstitutiveLawChecks > "${LOGFILE}" 2>&1 || status=$?

grep -E "^(PASS|FAIL):" "${LOGFILE}" || true

# Every law must have been checked, not merely every check that ran must have
# passed: a law that stopped the run, or was dropped from the list, would
# otherwise leave fewer lines and nothing failing
nLaws=$(grep -c "^Law checks: [A-Za-z]" "${LOGFILE}" || true)
nExpected=$(grep -E "^[A-Za-z]+$" constant/lawChecks | grep -vc "^FoamFile$" || true)
summary=$(grep -E "^Law checks: [0-9]+ passed" "${LOGFILE}" || true)

echo
if [[ ${status} -ne 0 || -z "${summary}" || ${nLaws} -ne ${nExpected} ]]; then
    echo "FAIL: ${summary:-the application did not finish} (${nLaws} of ${nExpected} laws checked)"
    tail -n 20 "${LOGFILE}"
    echo "============================================================"
    echo "Regression test FAILED"
    echo "============================================================"
    exit 1
fi

echo "PASS: ${summary}, ${nLaws} laws"
echo "============================================================"
echo "Regression test PASSED"
echo "============================================================"
