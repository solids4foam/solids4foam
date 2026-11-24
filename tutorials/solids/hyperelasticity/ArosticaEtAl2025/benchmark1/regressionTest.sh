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
# the mixed displacement-pressure formulation, to t = 0.25 s: through the
# passive filling phase and past the onset of the active stress, which
# switches on at t = 0.1645 s. This exercises the viscousHolzapfelOgden law
# (fibre and sheet families, viscosity, the volumetric split), the tabulated
# active tension of electroMechanicalLaw, the spring-dashpot conditions, the
# follower pressure and setFibreField's ellipsoid method.
#
# p0 x is checked at two times against reference values of this test, which
# are mesh-1 values and not the benchmark answer:
#   t = 0.05 s, filling, passive response only (participants: 0.370 to
#       0.504 mm);
#   t = 0.25 s, contraction, active stress on (participants: -28.1 to
#       -30.8 mm).
# A sign error, a missing 1/J or a missing active term moves the second by
# far more than its tolerance.
#
# Created by Claude Opus 5.5, supervised by Chanikya Valeti, UCD.
# ============================================================

END_TIME=0.25

# Reference p0 x [m] at the check times, and the relative tolerance
CHECK_TIMES=(0.05 0.25)
P0X_REF=(3.5960e-4 -3.3893e-2)
REL_TOL=0.02

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
DISP_FILE="postProcessing/0/solidPointDisplacement_p0.dat"

echo "============================================================"
echo "Arostica et al. (2025) benchmark 1 regression test"
echo "Step 1, mesh level 1, to t = ${END_TIME} s:"
for i in "${!CHECK_TIMES[@]}"; do
    echo "p0 x at t = ${CHECK_TIMES[$i]} s: ${P0X_REF[$i]} m, within ${REL_TOL}"
done
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
"${SOLIDS4FOAM_SED}" -i "s/^writeInterval .*/writeInterval   ${END_TIME};/" \
    "${CASE_DIR}/system/controlDict"

# The plot script and the data cache shared by the benchmarks: next to this
# case, or, when this case is a copy (as Alltest-regression makes), in the
# source tree, which is at the same depth below the repository root
ARO_DIR="${SCRIPT_DIR}/.."
if [[ ! -f "${ARO_DIR}/plotVsReference.gnuplot" ]]; then
    ARO_DIR="${SCRIPT_DIR}/../../../../../tutorials/solids/hyperelasticity/ArosticaEtAl2025"
fi
export ARO_DIR

( cd "${CASE_DIR}" && MESH_NUMBER=1 ./Allrun step1 pressure \
    > "${ALLRUN_LOGFILE}" 2>&1 ) || true

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped here"
    exit 0
fi

# Not a skip: reaching here without a displacement history means the case was
# expected to run and did not. The end of log.Allrun is printed, so that the
# cause is in the test output itself (e.g. in CI, where the logs may not be
# kept)
if [[ ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    if grep -qs "^Could not download" "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "FAIL: the mesh archive could not be downloaded"
    else
        echo "FAIL: no displacement history, and the tutorial did not report a skip"
    fi
    echo "      (see ${CASE_DIR}/${ALLRUN_LOGFILE}, which ends with:)"
    tail -n 15 "${CASE_DIR}/${ALLRUN_LOGFILE}" 2>/dev/null | sed 's/^/      | /'
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

if ! awk "BEGIN {exit !((${t} - ${END_TIME})^2 <= 1e-12)}"; then
    echo "FAIL: the displacement history stops at t = ${t}"
    failures=$((failures + 1))
else
    for i in "${!CHECK_TIMES[@]}"; do
        tc="${CHECK_TIMES[$i]}"
        ref="${P0X_REF[$i]}"
        p0x=$(awk -v tc="${tc}" \
            '!/^#/ && ($1 - tc)^2 <= 1e-12 {v = $2} END {print v}' \
            "${CASE_DIR}/${DISP_FILE}")

        if [[ -z "${p0x}" ]]; then
            echo "FAIL: no p0 displacement at t = ${tc} s"
            failures=$((failures + 1))
        elif awk -v v="${p0x}" -v r="${ref}" -v tol="${REL_TOL}" \
            'BEGIN {d = v - r; if (d < 0) d = -d; a = (r < 0 ? -r : r);
                    exit !(d <= tol*a)}'
        then
            printf "PASS: p0 x = %.6g m at t = %s s (reference %s m)\n" \
                "${p0x}" "${tc}" "${ref}"
        else
            printf "FAIL: p0 x = %.6g m at t = %s s (reference %s m)\n" \
                "${p0x}" "${tc}" "${ref}"
            failures=$((failures + 1))
        fi
    done
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
