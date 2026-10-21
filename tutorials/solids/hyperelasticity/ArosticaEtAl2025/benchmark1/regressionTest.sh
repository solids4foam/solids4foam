#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/step1"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# GNU sed, for the in-place edits below
solids4Foam::requireGnuSed

# ============================================================
# Arostica et al. (2025) benchmark 1 regression test
#
# Step 1 (active and passive response) on the coarsest mesh (level 1), with
# the mixed displacement-pressure formulation, for the first 50 time steps
# (t = 0.05 s, the start of the filling phase). This exercises the
# viscoActiveHolzapfelOgden law (fibre and sheet families, viscosity, the
# volumetric split), the spring-dashpot conditions, the follower pressure and
# setFibreField's ellipsoid method, in about a minute.
#
# The check is a band, not a reference value: p0 x at t = 0.05 s must lie
# near the benchmark participants' range, 0.370 to 0.504 mm (about 0.51 mm
# with the solidSpringDashpot condition of #528 on this mesh). It is still
# sensitive to the momentum stabilisation scale; a tight reference will be
# added once that is settled (stabilisation and mesh study, see README.md).
#
# Created by Claude Opus 5.5, supervised by Chanikya Valeti, UCD.
# ============================================================

END_TIME=0.05
P0X_MIN=3.0e-4
P0X_MAX=6.0e-4

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
DISP_FILE="postProcessing/0/solidPointDisplacement_p0.dat"

echo "============================================================"
echo "Arostica et al. (2025) benchmark 1 regression test"
echo "Step 1, mesh level 1, t = ${END_TIME} s:"
echo "p0 x in [${P0X_MIN}, ${P0X_MAX}] m"
echo "============================================================"
echo

rm -rf "${CASE_DIR}"
mkdir -p "${CASE_DIR}"
for item in "${SCRIPT_DIR}"/*; do
    [[ "$(basename "${item}")" == "regressionTests" ]] && continue
    cp -a "${item}" "${CASE_DIR}/"
done

"${SOLIDS4FOAM_SED}" -i "s/^endTime .*/endTime         ${END_TIME};/" \
    "${CASE_DIR}/system/controlDict"
"${SOLIDS4FOAM_SED}" -i 's/^writePrecision.*/writePrecision  14;/' \
    "${CASE_DIR}/system/controlDict"

# The Python scripts shared by the benchmarks: next to this case, or, when
# this case is a copy (as Alltest-regression makes), in the source tree, which
# is at the same depth below the repository root
ARO_SCRIPTS_DIR="${SCRIPT_DIR}/../scripts"
if [[ ! -d "${ARO_SCRIPTS_DIR}" ]]; then
    ARO_SCRIPTS_DIR="${SCRIPT_DIR}/../../../../../tutorials/solids/hyperelasticity/ArosticaEtAl2025/scripts"
fi
export ARO_SCRIPTS_DIR

( cd "${CASE_DIR}" && MESH_NUMBER=1 ./Allrun step1 pressure \
    > "${ALLRUN_LOGFILE}" 2>&1 ) || true

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped here"
    exit 0
fi

# Not a skip: reaching here without a displacement history means the case was
# expected to run and did not
if [[ ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    echo "FAIL: no displacement history, and the tutorial did not report a skip"
    echo "      (see ${CASE_DIR}/${ALLRUN_LOGFILE})"
    exit 1
fi

failures=0

# The fibre field must exist and be set, or the anisotropic law would be
# checked against an isotropic material
if [[ ! -f "${CASE_DIR}/0/f0" ]]; then
    echo "FAIL: setFibreField produced no fibre field"
    failures=$((failures + 1))
elif solids4Foam::internalFieldValues "${CASE_DIR}/0/f0" \
    | awk '{v = ($1 < 0 ? -$1 : $1); if (v > 1e-8) set = 1} END {exit !set}'
then
    echo "PASS: the fibre field is set"
else
    echo "FAIL: the fibre field is everywhere zero"
    failures=$((failures + 1))
fi

# The run must have finished, with every Newton solve converged
if ! grep -q "^End" "${CASE_DIR}/${SOLVER_LOGFILE}"; then
    echo "FAIL: the solver did not run to completion"
    failures=$((failures + 1))
elif grep -qE "Nonlinear solve did not converge|SNES convergence error" \
    "${CASE_DIR}/${SOLVER_LOGFILE}"
then
    echo "FAIL: a Newton solve did not converge"
    failures=$((failures + 1))
else
    echo "PASS: the solver ran to completion"
fi

t=$(awk 'END {print $1}' "${CASE_DIR}/${DISP_FILE}")
p0x=$(awk 'END {print $2}' "${CASE_DIR}/${DISP_FILE}")

if [[ -z "${p0x}" ]]; then
    echo "FAIL: could not read the p0 displacement"
    failures=$((failures + 1))
elif ! awk "BEGIN {exit !((${t} - ${END_TIME})^2 <= 1e-12)}"; then
    echo "FAIL: the displacement history stops at t = ${t}"
    failures=$((failures + 1))
elif awk "BEGIN {exit !(${p0x} >= ${P0X_MIN} && ${p0x} <= ${P0X_MAX})}"; then
    printf "PASS: p0 x = %.6g m\n" "${p0x}"
else
    printf "FAIL: p0 x = %.6g m\n" "${p0x}"
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
