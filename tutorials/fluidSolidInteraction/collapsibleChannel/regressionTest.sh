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

# Wall-midpoint vertical displacement at the first trough and at t = 1 s
REF_TROUGH=-0.21293
REF_END=-0.14786
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

if [[ ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    echo "FAIL: Could not find ${DISP_FILE}"
    exit 1
fi

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

trough=$(awk '!/^#/ && $1 >= 0.3 && $1 <= 0.75 {
    if (min == "" || $3 < min) min = $3
} END { print min }' "${CASE_DIR}/${DISP_FILE}")

end_disp=$(awk '!/^#/ { t = $1; v = $3 } END {
    if (t > '"${END_TIME}"' - 1e-6) print v
}' "${CASE_DIR}/${DISP_FILE}")

if [[ -z "${trough}" || -z "${end_disp}" ]]; then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

failures=0

check() {
    local label=$1 value=$2 reference=$3
    local diff
    diff=$(abs "$(awk "BEGIN {print ${value} - ${reference}}")")
    if awk "BEGIN {exit !(${diff} < ${DISP_TOL})}"; then
        printf "PASS: %s = %.6g (reference %.6g, diff %.3g)\n" \
            "${label}" "${value}" "${reference}" "${diff}"
    else
        printf "FAIL: %s = %.6g (reference %.6g, diff %.3g)\n" \
            "${label}" "${value}" "${reference}" "${diff}"
        failures=$((failures + 1))
    fi
}

check "trough displacement" "${trough}" "${REF_TROUGH}"
check "displacement at t = ${END_TIME} s" "${end_disp}" "${REF_END}"

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
