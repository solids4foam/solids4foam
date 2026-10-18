#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# channelLeaflet FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression tolerances and reference values
# ------------------------------------------------------------

# The case is run to t = 1 s, which covers the ramp of the inflow and the
# first large deflection of the leaflet
END_TIME=1
DELTA_T=0.01

# Tip displacement at t = 0.5 s and t = 1 s
REF_DX_HALF=0.20743
REF_DX_END=0.38140
REF_DY_END=-0.21541
DISP_TOL=0.002

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_tip.dat"

echo "============================================================"
echo "channelLeaflet FSI regression test"
echo "Tip displacement difference < ${DISP_TOL}"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    # Copy the inputs only: no results of a previous run of the tutorial
    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        case "${base_item}" in
            regressionTests|verification|postProcessing|dynamicCode|log.*)
                continue
                ;;
            processor*)
                continue
                ;;
        esac
        if [[ "${base_item}" != "0" \
           && "${base_item}" =~ ^[0-9.eE+-]+$ ]]; then
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
    -v refDxHalf="${REF_DX_HALF}" \
    -v refDxEnd="${REF_DX_END}" \
    -v refDyEnd="${REF_DY_END}" \
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
    if (NF < 3 || !finite($1) || !finite($2) || !finite($3)) {
        printf "FAIL: non-finite or truncated sample at line %d: %s\n", NR, $0
        bad = 1
        next
    }
    t = $1 + 0; dx = $2 + 0; dy = $3 + 0
    # Reject values that overflowed to infinity when converted
    if (absval(t) > 1e30 || absval(dx) > 1e30 || absval(dy) > 1e30) {
        printf "FAIL: non-finite sample at line %d: %s\n", NR, $0
        bad = 1
        next
    }
    if (t > 0) n++
    if (absval(t - 0.5) < 0.25*deltaT) { dxHalf = dx; halfSet = 1 }
    last = t; lastDx = dx; lastDy = dy
}
END {
    if (bad) exit 1
    expected = int(endTime/deltaT + 0.5)
    if (n < expected || halfSet == 0 || absval(last - endTime) > 1e-6) {
        printf "FAIL: incomplete history (%d of %d samples, last at t = %g)\n", n, expected, last
        exit 1
    }
    report("tip x displacement at t = 0.5 s", dxHalf, refDxHalf)
    report("tip x displacement at t = " endTime " s", lastDx, refDxEnd)
    report("tip y displacement at t = " endTime " s", lastDy, refDyEnd)
    exit failed ? 1 : 0
}' "${CASE_DIR}/${DISP_FILE}" || failures=1

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
