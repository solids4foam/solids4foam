#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# cylinderExpansion regression test
# Checks the final inner-boundary radial stress history.
# ============================================================

SIGMA_MIN=0.0115
SIGMA_MAX=0.0120

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "cylinderExpansion regression test"
echo "Final inner radial stress in [${SIGMA_MIN}, ${SIGMA_MAX}]"
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

value_file=$(find "${CASE_DIR}/postProcessing" -name 'solidStressesinner.dat' -print | tail -n 1)
if [[ -z "${value_file}" ]]; then
    echo "FAIL: Could not find inner stress history"
    exit 1
fi

sigma_rr=$(awk 'END {print -1e-6*$5}' "${value_file}")
if [[ -z "${sigma_rr}" ]]; then
    echo "FAIL: Could not extract final inner radial stress"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${sigma_rr} >= ${SIGMA_MIN} && ${sigma_rr} <= ${SIGMA_MAX})}"; then
    printf "PASS: Final inner radial stress = %.6g\n" "${sigma_rr}"
else
    printf "FAIL: Final inner radial stress = %.6g\n" "${sigma_rr}"
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
