#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# cylinderCrush regression test
# Uses the short displacement and force histories as a contact benchmark.
# ============================================================

FORCE_Y_MIN=-550
FORCE_Y_MAX=-545
DISP_Y_MIN=-0.0034
DISP_Y_MAX=-0.0032

ALLRUN_LOGFILE="log.Allrun"
FORCE_FILE="postProcessing/0/solidForcescylinderContact.dat"
DISP_FILE="postProcessing/0/solidPointDisplacement_displacement.dat"

echo "============================================================"
echo "cylinderCrush regression test"
echo "Final cylinder force_y in [${FORCE_Y_MIN}, ${FORCE_Y_MAX}] N"
echo "Final probe disp_y in [${DISP_Y_MIN}, ${DISP_Y_MAX}] m"
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

    sed -i.bak 's/^endTime[[:space:]]\+30;/endTime         1;/' "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"
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

if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" || ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    echo "FAIL: Could not find one or more history files"
    exit 1
fi

final_force_y=$(awk 'END {print $3}' "${CASE_DIR}/${FORCE_FILE}")
final_disp_y=$(awk 'END {print $3}' "${CASE_DIR}/${DISP_FILE}")

if [[ -z "${final_force_y}" || -z "${final_disp_y}" ]]; then
    echo "FAIL: Could not extract final force/displacement"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${final_force_y} >= ${FORCE_Y_MIN} && ${final_force_y} <= ${FORCE_Y_MAX})}"; then
    printf "PASS: Final force_y = %.6g\n" "${final_force_y}"
else
    printf "FAIL: Final force_y = %.6g\n" "${final_force_y}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_disp_y} >= ${DISP_Y_MIN} && ${final_disp_y} <= ${DISP_Y_MAX})}"; then
    printf "PASS: Final disp_y = %.6g\n" "${final_disp_y}"
else
    printf "FAIL: Final disp_y = %.6g\n" "${final_disp_y}"
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
