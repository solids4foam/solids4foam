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
# squarePlate regression test
# Checks the peak transverse deflection written to wVf.
# ============================================================

WF_MIN=6.8e-4
WF_MAX=7.0e-4

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "squarePlate regression test"
echo "Max wVf in [${WF_MIN}, ${WF_MAX}]"
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

field_file=$(find "${CASE_DIR}" -path '*/1/wVf' -print | tail -n 1)
if [[ -z "${field_file}" ]]; then
    echo "FAIL: Could not find wVf field output"
    exit 1
fi

max_wvf=$(awk '
    BEGIN {inlist=0; max=""}
    /^\($/ {inlist=1; next}
    /^\)$/ {inlist=0; next}
    inlist && $1 ~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/ {
        if (max == "" || $1 > max) {
            max = $1
        }
    }
    END {print max}
' "${field_file}")

if [[ -z "${max_wvf}" ]]; then
    echo "FAIL: Could not extract max wVf"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${max_wvf} >= ${WF_MIN} && ${max_wvf} <= ${WF_MAX})}"; then
    printf "PASS: Max wVf = %.6g\n" "${max_wvf}"
else
    printf "FAIL: Max wVf = %.6g\n" "${max_wvf}"
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
