#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# immersedSurfaceTraction quadrature checks
# Runs Test-immersedSurfaceTraction on a uniform and a graded
# two-dimensional mesh and a three-dimensional mesh (no flow solution): the
# quadratures none and clustered of the surface traction against the exact
# integrals of constant and linear pressures on closed surfaces, in three
# facet orders, static, rotated, translated and bent. See README.md
# ============================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LOGFILE="log.Test-immersedSurfaceTraction"

# Source solids4Foam scripts. From PATH, as a tutorial's Allrun does, since
# this case is run both in place and from a copy at a different depth
source solids4FoamScripts.sh

echo "============================================================"
echo "immersedSurfaceTraction quadrature checks"
echo "============================================================"

# The immersed boundary library is only built with OpenFOAM.com
solids4Foam::caseDoesNotRunWithFoamExtend
solids4Foam::caseDoesNotRunWithOpenFOAMOrg

# A skip where the application is not built, and a failure in CI, where it
# always is
solids4Foam::requireTestApp Test-immersedSurfaceTraction \
    || exit $(( $? - 1 ))

# Number of checks of each mesh: 4 two-dimensional surfaces x (3 rigid
# motions x (3 facet orders x 12 + 2) + 1 bend x (3 x 11 + 2)), on each
# two-dimensional mesh, and the three-dimensional ellipsoid x (2 rigid
# motions x (3 x 11 + 2) + 1 bend x (3 x 10 + 2)): a check that stopped, or a
# surface or motion dropped from the application, fails the test
declare -A EXPECTED=([twoD]=596 [twoDGraded]=596 [threeD]=102)

# On the graded mesh the least squares fits are solved with a larger rounding
# error: the tolerances of the invariants that are exact up to this rounding
# are multiplied by 1e3 (the largest errors are about 4e-9 of the 1e-7
# tolerance of the traction and 4e-11 of the 1e-9 tolerance of the forces)
declare -A OPTIONS=([twoD]="" [twoDGraded]="-roundoff 1e3" [threeD]="")

failed=0
for sub in twoD twoDGraded threeD
do
    CASE_DIR="${SCRIPT_DIR}/regressionTests/${sub}"
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}/constant"
    cp -a "${SCRIPT_DIR}/${sub}/system" "${CASE_DIR}/"
    (
        cd "${CASE_DIR}"
        blockMesh > log.blockMesh 2>&1
        status=0
        IFS=' ' read -r -a options <<< "${OPTIONS[${sub}]}"
        Test-immersedSurfaceTraction ${options[@]+"${options[@]}"} \
            > "${LOGFILE}" 2>&1 || status=$?
        summary=$(grep -E "^Checks: [0-9]+ passed" "${LOGFILE}" || true)
        nPass=$(grep -c "^PASS:" "${LOGFILE}" || true)
        grep -E "^FAIL:" "${LOGFILE}" || true
        if [[ ${status} -ne 0 || -z "${summary}" \
           || ${nPass} -ne ${EXPECTED[${sub}]} ]]
        then
            echo "FAIL: ${sub}: ${summary:-the application did not finish}" \
                "(${nPass} of ${EXPECTED[${sub}]} checks passed)"
            tail -n 20 "${LOGFILE}"
            exit 1
        fi
        echo "PASS: ${sub}: ${summary}"
    ) || failed=1
done

echo "============================================================"
if [[ ${failed} -ne 0 ]]
then
    echo "Regression test FAILED"
    echo "============================================================"
    exit 1
fi
echo "Regression test PASSED"
echo "============================================================"
