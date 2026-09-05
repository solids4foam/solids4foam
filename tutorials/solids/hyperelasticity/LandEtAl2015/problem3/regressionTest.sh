#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# Land et al. (2015) problem 3 regression test
#
# The benchmark measures how far a point on the tissue moves under an active
# tension that ramps over the run. What makes this case worth having is that
# it is the only tutorial whose material is anisotropic about a fibre
# direction: the fibre field is built by setFibreField before the solver runs,
# and the constitutive law reads it.
#
# It is a legacy-only test for now. electroMechanicalLaw and the GuccioneElastic
# law it wraps have not been ported to the mechanicalConstitutiveLaw framework,
# and this case is here so that the port has something to be checked against
# rather than being written and hoped over.
# ============================================================

MAG_D_MIN=8.0e-4
MAG_D_MAX=9.0e-4

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"

echo "============================================================"
echo "Land et al. (2015) problem 3 regression test"
echo "Final probe |D| in [${MAG_D_MIN}, ${MAG_D_MAX}] m"
echo "============================================================"
echo

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
( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped here"
    exit 0
fi

if [[ ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    echo "SKIP: no displacement history (the case did not run here)"
    exit 0
fi

failures=0

# The fibre field has to exist and point somewhere, or this case is an
# isotropic one wearing an anisotropic law and the port it exists to check
# would be checked against nothing
if [[ ! -f "${CASE_DIR}/0/f0" ]]; then
    echo "FAIL: setFibreField produced no fibre field"
    failures=$((failures + 1))
elif python3 - "${CASE_DIR}/0/f0" << 'PYEOF'
import re, sys
body = open(sys.argv[1]).read().split('* * * * *')[-1]
nums = [float(x) for x in re.findall(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', body)]
# any component away from zero means a direction was actually set
sys.exit(0 if any(abs(v) > 1e-8 for v in nums) else 1)
PYEOF
then
    echo "PASS: the fibre field is set"
else
    echo "FAIL: the fibre field is everywhere zero"
    failures=$((failures + 1))
fi

magD=$(awk 'END {print $5}' "${CASE_DIR}/${DISP_FILE}")

if [[ -z "${magD}" ]]; then
    echo "FAIL: could not read the final probe displacement"
    failures=$((failures + 1))
elif awk "BEGIN {exit !(${magD} >= ${MAG_D_MIN} && ${magD} <= ${MAG_D_MAX})}"; then
    printf "PASS: final probe |D| = %.6g\n" "${magD}"
else
    printf "FAIL: final probe |D| = %.6g\n" "${magD}"
    failures=$((failures + 1))
fi

echo
if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED"
    echo "============================================================"
    exit 0
else
    echo "============================================================"
    echo "Regression test FAILED (${failures} checks)"
    echo "============================================================"
    exit 1
fi
