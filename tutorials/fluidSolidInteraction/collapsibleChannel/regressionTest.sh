#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# collapsibleChannel FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression tolerances and reference values
# ------------------------------------------------------------

# The case is run to t = 1 s, which covers the collapse of the wall, its first
# trough and the following rebound
END_TIME=1
DELTA_T=0.025

# Wall-midpoint vertical displacement at the first trough and at t = 1 s
REF_TROUGH=-0.21293
REF_END=-0.14756
DISP_TOL=0.002

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_wallMid.dat"

echo "============================================================"
echo "collapsibleChannel FSI regression test"
echo "Wall-midpoint displacement difference < ${DISP_TOL}"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" \
           || "${base_item}" == "verification" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done

    sed -i.bak "s/^endTime .*/endTime         ${END_TIME};/" \
        "${CASE_DIR}/system/controlDict"
}

# ------------------------------------------------------------
# Clean & run case
# ------------------------------------------------------------

CHECK_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run)
            CHECK_ONLY=true
            ;;
        *)
            ;;
    esac
done

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# The case requires PETSc: without it Allrun writes a declared skip message,
# which is the only valid reason to skip the checks
if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

if ! grep -q "^End[[:space:]]*$" "${CASE_DIR}/log.solids4Foam" 2>/dev/null; then
    echo "FAIL: solids4Foam did not finish (see ${CASE_DIR}/log.solids4Foam)"
    exit 1
fi

if [[ ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    echo "FAIL: Could not find ${DISP_FILE}"
    exit 1
fi

# ------------------------------------------------------------
# Extract and check values
# ------------------------------------------------------------

# Every sample must be a finite number, and the history must reach END_TIME.
# All values are passed to awk as data, never as program text
failures=0
awk \
    -v endTime="${END_TIME}" \
    -v deltaT="${DELTA_T}" \
    -v refTrough="${REF_TROUGH}" \
    -v refEnd="${REF_END}" \
    -v tol="${DISP_TOL}" '
function finite(x) {
    return x ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/
}
function absval(x) { return x < 0 ? -x : x }
function report(label, value, reference,    diff) {
    diff = absval(value - reference)
    if (diff < tol) {
        printf "PASS: %s = %.6g (reference %.6g, diff %.3g)\n", label, value, reference, diff
    } else {
        printf "FAIL: %s = %.6g (reference %.6g, diff %.3g)\n", label, value, reference, diff
        failed++
    }
}
/^#/ { next }
{
    if (!finite($1) || !finite($3)) {
        printf "FAIL: non-finite sample at line %d: %s\n", NR, $0
        bad = 1
        next
    }
    t = $1 + 0; v = $3 + 0
    # Reject values that overflowed to infinity when converted
    if (absval(t) > 1e30 || absval(v) > 1e30) {
        printf "FAIL: non-finite sample at line %d: %s\n", NR, $0
        bad = 1
        next
    }
    if (t > 0) n++
    if (t >= 0.3 && t <= 0.75 && (troughSet == 0 || v < trough)) {
        trough = v; troughSet = 1
    }
    last = t; lastValue = v
}
END {
    if (bad) exit 1
    expected = int(endTime/deltaT + 0.5)
    if (n < expected || troughSet == 0 || absval(last - endTime) > 1e-6) {
        printf "FAIL: incomplete history (%d of %d samples, last at t = %g)\n", n, expected, last
        exit 1
    }
    report("trough displacement", trough, refTrough)
    report("displacement at t = " endTime " s", lastValue, refEnd)
    exit failed ? 1 : 0
}' "${CASE_DIR}/${DISP_FILE}" || failures=1

# ------------------------------------------------------------
# Quasi-monolithic coupling with the coloured finite-difference solid block
# ------------------------------------------------------------
# A short run (t = 0.15 s, 24 steps of 0.00625 s) of the quasi-monolithic
# coupling with "solidJacobian colouredFiniteDifference": the wall-midpoint
# displacement must match the reference (the converged solution does not
# depend on the solid block), and the mean number of Krylov iterations per
# Newton iteration of the coupled system must stay small (about 20 with the
# finite-difference block, more than 300 with the compact block)
QM_CASE_DIR="${REGRESSION_ROOT}/monolithicColouredFD"
QM_END_TIME=0.15
QM_DELTA_T=0.00625
QM_REF_DISP=-0.0084436
QM_DISP_TOL=2e-6
QM_MAX_KRYLOV=60

if [ "$CHECK_ONLY" = false ]; then
    rm -rf "${QM_CASE_DIR}"
    mkdir -p "${QM_CASE_DIR}"
    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" \
           || "${base_item}" == "verification" ]]; then
            continue
        fi
        cp -a "${item}" "${QM_CASE_DIR}/"
    done
    sed -i.bak \
        -e "s/^endTime .*/endTime         ${QM_END_TIME};/" \
        -e "s/^deltaT .*/deltaT          ${QM_DELTA_T};/" \
        -e "s/^writeInterval .*/writeInterval   ${QM_END_TIME};/" \
        "${QM_CASE_DIR}/system/controlDict"
    sed -i.bak \
        "s/^\([[:space:]]*coupled[[:space:]]*yes;\)/\1\n    solidJacobian colouredFiniteDifference;/" \
        "${QM_CASE_DIR}/constant/fsiProperties.monolithic"
    ( cd "${QM_CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${QM_CASE_DIR}" && ./Allrun monolithic > "${ALLRUN_LOGFILE}" 2>&1 ) || true
fi

if solids4Foam::regressionCaseSkipped "${QM_CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping the quasi-monolithic checks: the tutorial skipped in this environment"
elif ! grep -q "^End[[:space:]]*$" "${QM_CASE_DIR}/log.solids4Foam" 2>/dev/null; then
    echo "FAIL: quasi-monolithic solids4Foam did not finish (see ${QM_CASE_DIR}/log.solids4Foam)"
    failures=1
else
    awk \
        -v endTime="${QM_END_TIME}" \
        -v ref="${QM_REF_DISP}" \
        -v tol="${QM_DISP_TOL}" '
    /^#/ { next }
    { t = $1; v = $3 }
    END {
        d = v - ref; if (d < 0) d = -d
        if ((t - endTime)^2 > 1e-12 || d >= tol) {
            printf "FAIL: QM wall displacement at t = %g s = %.7g (reference %.7g)\n", t, v, ref
            exit 1
        }
        printf "PASS: QM wall displacement at t = %g s = %.7g (reference %.7g, diff %.3g)\n", t, v, ref, d
    }' "${QM_CASE_DIR}/${DISP_FILE}" || failures=1

    awk -v maxK="${QM_MAX_KRYLOV}" '
    /Solving the monolithic momentum-continuity system/ { coupled = 1; next }
    /Solving the momentum equation for D/ { coupled = 0 }
    coupled && /Linear solve (converged|did not converge) due to/ {
        sum += $NF; n++
    }
    END {
        if (n == 0) { print "FAIL: no coupled linear solves found"; exit 1 }
        mean = sum/n
        if (mean > maxK) {
            printf "FAIL: mean Krylov iterations per Newton iteration %.1f > %g\n", mean, maxK
            exit 1
        }
        printf "PASS: mean Krylov iterations per Newton iteration %.1f (limit %g)\n", mean, maxK
    }' "${QM_CASE_DIR}/log.solids4Foam" || failures=1
fi

echo
if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED"
    echo "============================================================"
    exit 0
else
    echo "============================================================"
    echo "Regression test FAILED"
    echo "============================================================"
    exit 1
fi
