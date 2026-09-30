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
# staticCylinderInChannel regression test
# Runs the coarsest mesh (MESH_LEVEL=1) to steady state and
# checks the drag and lift coefficients of the immersed
# cylinder against the values of this method on this mesh.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

REG_END_TIME=5
CD_MIN=6.00
CD_MAX=6.08
CL_MIN=0.010
CL_MAX=0.020

ALLRUN_LOGFILE="log.Allrun"
FORCE_FILE="postProcessing/immersedBoundary/0/cylinder.dat"

echo "============================================================"
echo "staticCylinderInChannel regression test"
echo "Regression end time = ${REG_END_TIME}"
echo "Final Cd in [${CD_MIN}, ${CD_MAX}]"
echo "Final Cl in [${CL_MIN}, ${CL_MAX}]"
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

    sed -i.bak \
        "s/^endTime[[:space:]]\+[0-9.]*;/endTime         ${REG_END_TIME};/" \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"
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
    ( cd "${CASE_DIR}" && MESH_LEVEL=1 ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
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

if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
    echo "FAIL: Could not find ${FORCE_FILE}"
    exit 1
fi

# Cd = 2 Fx/(rho Umean^2 D Lz) = 5000 Fx, and likewise for Cl
final_cd=$(awk '!/^#/ { cd = 5000*$2 } END { print cd }' "${CASE_DIR}/${FORCE_FILE}")
final_cl=$(awk '!/^#/ { cl = 5000*$3 } END { print cl }' "${CASE_DIR}/${FORCE_FILE}")

if [[ -z "${final_cd}" || -z "${final_cl}" ]]; then
    echo "FAIL: Could not extract the final drag and lift coefficients"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check_range() {
    local name=$1 value=$2 min=$3 max=$4
    if awk "BEGIN {exit !(${value} >= ${min} && ${value} <= ${max})}"; then
        printf "PASS: %s = %.6g\n" "${name}" "${value}"
    else
        printf "FAIL: %s = %.6g\n" "${name}" "${value}"
        failures=$((failures + 1))
    fi
}

check_range "Final Cd" "${final_cd}" "${CD_MIN}" "${CD_MAX}"
check_range "Final Cl" "${final_cl}" "${CL_MIN}" "${CL_MAX}"

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
