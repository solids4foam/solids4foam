#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# circulationNetwork checks
# Runs Test-circulationNetworkChecks, which checks the lumped circulation
# network on small circuits with known answers: a port or a valve feeding an
# inductor, the waveform derivatives, and a hard valve opening from rest
# ============================================================

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
LOGFILE="log.Test-circulationNetworkChecks"

# Source solids4Foam scripts. From PATH, as a tutorial's Allrun does, since
# this case is run both in place and from a copy at a different depth
source solids4FoamScripts.sh

echo "============================================================"
echo "circulationNetwork checks"
echo "============================================================"

# The circulation models are only built for OpenFOAM.com
solids4Foam::caseDoesNotRunWithFoamExtend
solids4Foam::caseDoesNotRunWithOpenFOAMOrg

# A skip where the application is not built, and a failure in CI, where it
# always is
solids4Foam::requireTestApp Test-circulationNetworkChecks \
    || exit $(( $? - 1 ))

# Run in a copy: the application writes the waveform table it reads
CASE_DIR="${SCRIPT_DIR}/regressionTests/main"
rm -rf "${CASE_DIR}"
mkdir -p "${CASE_DIR}"
cd "${CASE_DIR}"

status=0
Test-circulationNetworkChecks > "${LOGFILE}" 2>&1 || status=$?

grep -E "^(PASS|FAIL):" "${LOGFILE}" || true

summary=$(grep -E "^Circulation network checks: [0-9]+ passed" "${LOGFILE}" || true)

echo
if [[ ${status} -ne 0 || -z "${summary}" ]]; then
    echo "FAIL: ${summary:-the application did not finish}"
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
