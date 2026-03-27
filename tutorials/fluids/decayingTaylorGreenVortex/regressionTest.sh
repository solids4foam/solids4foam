#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# decayingTaylorGreenVortex regression test
# Checks the final velocity error norms written by the case.
# ============================================================

L1_MAX=3e-4
L2_MAX=3.5e-4
LINF_MAX=1.1e-3

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "decayingTaylorGreenVortex regression test"
echo "Velocity error norms must satisfy:"
echo "  L1 <= ${L1_MAX}"
echo "  L2 <= ${L2_MAX}"
echo "  LInf <= ${LINF_MAX}"
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

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

log_file="${CASE_DIR}/log.solids4Foam"
if [[ ! -f "${log_file}" ]]; then
    echo "FAIL: Could not find solver log"
    exit 1
fi

norms=$(awk '
    /Velocity errors norms:/ {found=1; next}
    found && /mean L1 =/ {l1=$4; next}
    found && /mean L2 =/ {l2=$4; next}
    found && /LInf =/ {linf=$3; found=0}
    END {print l1, l2, linf}
' "${log_file}")

IFS=' ' read -r l1 l2 linf <<< "${norms}"

if [[ -z "${l1:-}" || -z "${l2:-}" || -z "${linf:-}" ]]; then
    echo "FAIL: Could not extract velocity error norms"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${l1} <= ${L1_MAX})}"; then
    printf "PASS: L1 = %.6g\n" "${l1}"
else
    printf "FAIL: L1 = %.6g\n" "${l1}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${l2} <= ${L2_MAX})}"; then
    printf "PASS: L2 = %.6g\n" "${l2}"
else
    printf "FAIL: L2 = %.6g\n" "${l2}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${linf} <= ${LINF_MAX})}"; then
    printf "PASS: LInf = %.6g\n" "${linf}"
else
    printf "FAIL: LInf = %.6g\n" "${linf}"
    failures=$((failures + 1))
fi

if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
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
