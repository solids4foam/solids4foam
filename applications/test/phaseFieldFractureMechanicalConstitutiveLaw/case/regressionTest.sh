#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Phase field fracture material-point checks
# Runs Test-phaseFieldFractureMechanicalConstitutiveLaw on a one-cell mesh,
# which drives linearElasticMisesPlasticPhaseFieldFracture over multi-step
# strain paths and checks the degraded stress, the crack driving force, the
# homogeneous brittle response against its analytical solution, and the law
# under logStrain
# ============================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LOGFILE="log.Test-phaseFieldFractureMechanicalConstitutiveLaw"

# Source solids4Foam scripts. From PATH, as a tutorial's Allrun does, since
# this case is run both in place and from a copy at a different depth
source solids4FoamScripts.sh

echo "============================================================"
echo "Phase field fracture material-point checks"
echo "============================================================"

# A skip where the application is not built, and a failure in CI, where it
# always is
solids4Foam::requireTestApp Test-phaseFieldFractureMechanicalConstitutiveLaw \
    || exit $(( $? - 1 ))

# Run on a copy, as the tutorials' regression tests do: converting the case
# format moves files, which must not happen to the case itself
CASE_DIR="${SCRIPT_DIR}/regressionTests/main"
rm -rf "${CASE_DIR}"
mkdir -p "${CASE_DIR}/0"
cp -a "${SCRIPT_DIR}/constant" "${SCRIPT_DIR}/system" "${CASE_DIR}/"
cd "${CASE_DIR}"

solids4Foam::convertCaseFormat . > /dev/null
blockMesh > log.blockMesh 2>&1

status=0
Test-phaseFieldFractureMechanicalConstitutiveLaw > "${LOGFILE}" 2>&1 || status=$?

grep -E "^(PASS|FAIL):" "${LOGFILE}" || true

summary=$(grep -E "^Phase field checks: [0-9]+ passed, 0 failed" "${LOGFILE}" \
    || true)

echo
if [[ ${status} -ne 0 || -z "${summary}" ]]; then
    echo "FAIL: the application failed or did not finish"
    tail -n 20 "${LOGFILE}"
    echo "============================================================"
    echo "Regression test FAILED"
    echo "============================================================"
    exit 1
fi

echo "PASS: ${summary}"
echo "============================================================"
echo "Regression test PASSED"
echo "============================================================"
