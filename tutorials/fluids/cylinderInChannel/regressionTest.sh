#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# cylinderInChannel regression test
# Uses the force history as a cheap fluid benchmark check.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

REG_END_TIME=5
DRAG_MIN=0.622
DRAG_MAX=0.624
LIFT_MIN=0.0105
LIFT_MAX=0.0115

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "cylinderInChannel regression test"
echo "Regression end time = ${REG_END_TIME}"
echo "Final drag in [${DRAG_MIN}, ${DRAG_MAX}]"
echo "Final lift in [${LIFT_MIN}, ${LIFT_MAX}]"
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

    sed -i.bak 's/^endTime[[:space:]]\+50;/endTime         5;/' "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"
}

find_force_file() {
    local candidate
    for candidate in \
        "${CASE_DIR}/postProcessing/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/forces/0/forces.dat" \
        "${CASE_DIR}/postProcessing/fluid/forces/0/force.dat" \
        "${CASE_DIR}/postProcessing/fluid/forces/0/forces.dat"
    do
        if [[ -f "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done
    return 1
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

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

force_file=""
if ! force_file=$(find_force_file); then
    echo "FAIL: Could not find force history output"
    exit 1
fi

final_drag=$(awk '
    END {
        gsub(/[()]/, "", $0)
        print $2 + $5
    }
' "${force_file}")
final_lift=$(awk '
    END {
        gsub(/[()]/, "", $0)
        print $3 + $6
    }
' "${force_file}")

if [[ -z "${final_drag}" || -z "${final_lift}" ]]; then
    echo "FAIL: Could not extract final drag/lift"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${final_drag} >= ${DRAG_MIN} && ${final_drag} <= ${DRAG_MAX})}"; then
    printf "PASS: Final drag = %.6g\n" "${final_drag}"
else
    printf "FAIL: Final drag = %.6g\n" "${final_drag}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_lift} >= ${LIFT_MIN} && ${final_lift} <= ${LIFT_MAX})}"; then
    printf "PASS: Final lift = %.6g\n" "${final_lift}"
else
    printf "FAIL: Final lift = %.6g\n" "${final_lift}"
    failures=$((failures + 1))
fi

# Clean case again
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
