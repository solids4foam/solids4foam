#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# neckingBar regression test
# Checks the final loading force after the necking curve.
# ============================================================

FORCE_MIN=206.8
FORCE_MAX=207.0

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "neckingBar regression test"
echo "Final loading force in [${FORCE_MIN}, ${FORCE_MAX}]"
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
    sed -i.bak 's/^endTime         1;/endTime         0.184;/' "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

force_file=$(find "${CASE_DIR}/postProcessing" -name 'solidForcesloading.dat' -print | tail -n 1)
if [[ -z "${force_file}" ]]; then
    echo "FAIL: Could not find loading force history"
    exit 1
fi

final_force=$(awk 'END {print $2}' "${force_file}")
if [[ -z "${final_force}" ]]; then
    echo "FAIL: Could not extract final loading force"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${final_force} >= ${FORCE_MIN} && ${final_force} <= ${FORCE_MAX})}"; then
    printf "PASS: Final loading force = %.6g\n" "${final_force}"
else
    printf "FAIL: Final loading force = %.6g\n" "${final_force}"
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
