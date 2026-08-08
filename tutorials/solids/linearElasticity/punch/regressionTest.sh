#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# punch regression test
# Checks the final Z displacement.
# ============================================================

PUNCH_DISP_Z_MIN=-0.00025
PUNCH_DISP_Z_MAX=-0.0002

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "punch regression test"
echo "Punch loading dispZ in [${PUNCH_DISP_Z_MIN}, ${PUNCH_DISP_Z_MAX}]"
echo "============================================================"
echo

if ! command -v cartesianMesh >/dev/null 2>&1; then
    echo "Skipping this case as cartesianMesh is not installed"
    exit 0
fi

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

punch_file=$(find "${CASE_DIR}/postProcessing" -name 'solidForcesDisplacementspunchLoading.dat' -print | tail -n 1)

if [[ -z "${punch_file}" ]]; then
    echo "FAIL: Could not find punch force/displacement output"
    exit 1
fi

punch_disp_z=$(awk 'END {print $4}' "${punch_file}")

if [[ -z "${punch_disp_z:-}" ]]; then
    echo "FAIL: Could not extract punch Z displacement"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${punch_disp_z} >= ${PUNCH_DISP_Z_MIN} && ${punch_disp_z} <= ${PUNCH_DISP_Z_MAX})}"; then
    printf "PASS: punchLoading dispZ = %.6g\n" "${punch_disp_z}"
else
    printf "FAIL: punchLoading dispZ = %.6g\n" "${punch_disp_z}"
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
