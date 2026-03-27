#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# curvedCantilever regression test
# Checks the sampled numerical stress against the analytical
# stress field written by the built-in function object.
# ============================================================

MAX_STRESS_ERROR=2000

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "curvedCantilever regression test"
echo "Max stress component difference <= ${MAX_STRESS_ERROR}"
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

sample_file=$(find "${CASE_DIR}/postProcessing" -name 'line_sigma_analyticalStress.xy' -print | tail -n 1)
if [[ -z "${sample_file}" ]]; then
    echo "FAIL: Could not find sampled stress output"
    exit 1
fi

max_error=$(awk '
{
    d = 0;
    for (i = 4; i <= 9; i++) {
        e = $(i) - $(i + 6);
        if (e < 0) {
            e = -e;
        }
        if (e > d) {
            d = e;
        }
    }
    if (d > m) {
        m = d;
    }
}
END {
    print m;
}' "${sample_file}")

if [[ -z "${max_error}" ]]; then
    echo "FAIL: Could not extract max stress error"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${max_error} <= ${MAX_STRESS_ERROR})}"; then
    printf "PASS: Max stress component difference = %.6g\n" "${max_error}"
else
    printf "FAIL: Max stress component difference = %.6g\n" "${max_error}"
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
