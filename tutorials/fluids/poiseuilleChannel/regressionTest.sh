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
# poiseuilleChannel regression test
# Checks the steady laminar channel flow driven by the
# meanVelocityForce fvOption against the analytical solution:
# pressure gradient 12*nu*Ubar/h^2 = 1.2 and maximum velocity
# 1.5*Ubar = 1.5
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

GRADP_MIN=1.194
GRADP_MAX=1.206
UBAR_MIN=0.9999
UBAR_MAX=1.0001
UMAX_MIN=1.494
UMAX_MAX=1.506

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"

echo "============================================================"
echo "poiseuilleChannel regression test"
echo "Final pressure gradient in [${GRADP_MIN}, ${GRADP_MAX}]"
echo "Final mean velocity in [${UBAR_MIN}, ${UBAR_MAX}]"
echo "Final maximum velocity in [${UMAX_MIN}, ${UMAX_MAX}]"
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

if [[ ! -f "${CASE_DIR}/${SOLVER_LOGFILE}" ]]; then
    echo "FAIL: Could not find ${SOLVER_LOGFILE}"
    exit 1
fi

# The meanVelocityForce fvOption reports, on each time step:
# "Pressure gradient source: uncorrected Ubar = <Ubar>, pressure gradient = <gradP>"
final_ubar=$(awk -F'[=,]' '
    /Pressure gradient source/ { ubar = $2 }
    END { gsub(/ /, "", ubar); print ubar }
' "${CASE_DIR}/${SOLVER_LOGFILE}")
final_gradp=$(awk -F'=' '
    /Pressure gradient source/ { gradp = $3 }
    END { gsub(/ /, "", gradp); print gradp }
' "${CASE_DIR}/${SOLVER_LOGFILE}")

latest_time=$(solids4Foam::latestTime "${CASE_DIR}")
final_umax=""
if [[ -n "${latest_time}" && -f "${CASE_DIR}/${latest_time}/U" ]]; then
    final_umax=$(awk '
        /^internalField/ { inField = 1; next }
        inField && /^\)/ { exit }
        inField && /^\(/ {
            gsub(/[()]/, "", $0)
            if (n == 0 || $1 > umax) { umax = $1 }
            n++
        }
        END { if (n > 0) { print umax } }
    ' "${CASE_DIR}/${latest_time}/U")
fi

if [[ -z "${final_ubar}" || -z "${final_gradp}" || -z "${final_umax}" ]]; then
    echo "FAIL: Could not extract the final pressure gradient and velocities"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check_range() {
    local name=$1 value=$2 min=$3 max=$4
    if awk "BEGIN {exit !(${value} >= ${min} && ${value} <= ${max})}"; then
        printf "PASS: %s = %.8g\n" "${name}" "${value}"
    else
        printf "FAIL: %s = %.8g\n" "${name}" "${value}"
        failures=$((failures + 1))
    fi
}

check_range "Final pressure gradient" "${final_gradp}" "${GRADP_MIN}" "${GRADP_MAX}"
check_range "Final mean velocity" "${final_ubar}" "${UBAR_MIN}" "${UBAR_MAX}"
check_range "Final maximum velocity" "${final_umax}" "${UMAX_MIN}" "${UMAX_MAX}"

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
