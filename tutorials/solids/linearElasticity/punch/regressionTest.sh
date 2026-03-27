#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# punch regression test
# Checks the final displacement and reaction force histories.
# ============================================================

PUNCH_DISP_Z_MIN=-1.7e-4
PUNCH_DISP_Z_MAX=-1.4e-4
SUPPORT_FORCE_Z_MIN=1.95e5
SUPPORT_FORCE_Z_MAX=2.05e5

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "punch regression test"
echo "Punch loading dispZ in [${PUNCH_DISP_Z_MIN}, ${PUNCH_DISP_Z_MAX}]"
echo "Support forceZ in [${SUPPORT_FORCE_Z_MIN}, ${SUPPORT_FORCE_Z_MAX}]"
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
support_file=$(find "${CASE_DIR}/postProcessing" -name 'solidForcesDisplacementscylinderFixed.dat' -print | tail -n 1)

if [[ -z "${punch_file}" || -z "${support_file}" ]]; then
    echo "FAIL: Could not find punch force/displacement outputs"
    exit 1
fi

punch_disp_z=$(awk 'END {print $4}' "${punch_file}")
support_force_z=$(awk 'END {print $7}' "${support_file}")

if [[ -z "${punch_disp_z:-}" || -z "${support_force_z:-}" ]]; then
    echo "FAIL: Could not extract punch regression quantities"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${punch_disp_z} >= ${PUNCH_DISP_Z_MIN} && ${punch_disp_z} <= ${PUNCH_DISP_Z_MAX})}"; then
    printf "PASS: punchLoading dispZ = %.6g\n" "${punch_disp_z}"
else
    printf "FAIL: punchLoading dispZ = %.6g\n" "${punch_disp_z}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${support_force_z} >= ${SUPPORT_FORCE_Z_MIN} && ${support_force_z} <= ${SUPPORT_FORCE_Z_MAX})}"; then
    printf "PASS: cylinderFixed forceZ = %.6g\n" "${support_force_z}"
else
    printf "FAIL: cylinderFixed forceZ = %.6g\n" "${support_force_z}"
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
