#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# curvedBeams regression test
# Uses the reaction force history as a cheap contact benchmark check.
# ============================================================

FORCE_Y_MIN=-17.8
FORCE_Y_MAX=-17.6

ALLRUN_LOGFILE="log.Allrun"
FORCE_FILE="postProcessing/0/solidForcesdisplacement.dat"

echo "============================================================"
echo "curvedBeams regression test"
echo "Final displacement-patch force_y in [${FORCE_Y_MIN}, ${FORCE_Y_MAX}] N"
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

    sed -i.bak 's/^endTime[[:space:]]\+31.5;/endTime         9;/' "${CASE_DIR}/system/controlDict"
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

if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
    echo "FAIL: Could not find ${FORCE_FILE}"
    exit 1
fi

final_force_y=$(awk 'END {print $3}' "${CASE_DIR}/${FORCE_FILE}")

if [[ -z "${final_force_y}" ]]; then
    echo "FAIL: Could not extract final force_y"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${final_force_y} >= ${FORCE_Y_MIN} && ${final_force_y} <= ${FORCE_Y_MAX})}"; then
    printf "PASS: Final force_y = %.6g\n" "${final_force_y}"
else
    printf "FAIL: Final force_y = %.6g\n" "${final_force_y}"
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
