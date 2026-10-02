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
# slabCooling regression test
# Unconstrained thermal contraction
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

# Stress should be ~0 (numerical noise only)
SIGMA_MAX=1e3      # Pa

# Strain should be O(1e-8)
EPS_MIN=1e-9
EPS_MAX=1e-7

# ------------------------------------------------------------
# Log files
# ------------------------------------------------------------

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

variant="openfoamcom"
if [[ -n "${FOAMEXTEND:-}" || "${WM_PROJECT_VERSION:-}" == "4.1" ]]; then
    variant="foamextend"
elif [[ "${WM_PROJECT_VERSION:-}" != *"v"* ]]; then
    variant="openfoamorg"
fi

if [[ "${variant}" != "openfoamcom" ]]; then
    SIGMA_MAX=1.05e3
fi

# The final D, as the max and mean component magnitude of its internal values
# as written. thermoMechanicalLaw is a composite: it owns a sub-law, delegates
# to it, and subtracts the thermal term, and the framework reproduced the
# removed legacy mechanicalModel's D field exactly, in every figure written, on
# every fork (mcl-stage8-coverage, c3a92b3d). These are OpenFOAM.com v2512's,
# the same as OpenFOAM.org 9's; foam-extend 4.1 gives a max 1.8e-3 smaller,
# and the tolerance, 3e-3 of the largest value, covers that
REF_D_MAX=0.0484734
REF_D_MEAN=0.0120319022347934
REF_D_REL_TOL=3e-3

echo "============================================================"
echo "slabCooling regression test"
echo "Max sigmaEq < ${SIGMA_MAX} Pa"
echo "epsilonEq order: ${EPS_MIN} < eps < ${EPS_MAX}"
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
# Clean & run
# ------------------------------------------------------------

prepare_case
( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )

check_against_reference() {
    if ! grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the case constructed no mechanical constitutive law"
        return 1
    fi

    local t end_time
    t=$(solids4Foam::latestTime "${CASE_DIR}")
    end_time=$(sed -n 's/^endTime[[:space:]]*\([^;]*\);.*/\1/p' \
        "${CASE_DIR}/system/controlDict")

    if [[ -z "${t}" || -z "${end_time}" ]] \
        || ! awk "BEGIN {exit !((${t} - ${end_time})^2 <= 1e-20)}"
    then
        echo "FAIL: the case stopped at '${t}', not at the end time '${end_time}'"
        return 1
    fi

    solids4Foam::checkFieldNorms "D" "${CASE_DIR}/${t}/D" \
        "${REF_D_MAX}" "${REF_D_MEAN}" "${REF_D_REL_TOL}"
}

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk -F '=' '{print $2}' \
        | tr -d '[:space:]'
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk -F '=' '{print $2}' \
        | tr -d '[:space:]'
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

epsilon=$(extract_max_epsilon)
sigma=$(extract_max_sigma)

if [[ -z "${epsilon}" || -z "${sigma}" ]]
then
    echo "FAIL: Could not extract epsilonEq or sigmaEq"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

# --- Stress check ------------------------------------------------------------

if awk "BEGIN {exit !(${sigma} < ${SIGMA_MAX})}"
then
    printf "PASS: Max sigmaEq = %.6g Pa\n" "${sigma}"
else
    printf "FAIL: Max sigmaEq = %.6g Pa\n" "${sigma}"
    failures=$((failures + 1))
fi

# --- Strain order-of-magnitude check ----------------------------------------

if awk "BEGIN {exit !(${epsilon} > ${EPS_MIN} && ${epsilon} < ${EPS_MAX})}"
then
    printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
else
    printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
    failures=$((failures + 1))
fi

if ! check_against_reference; then
    failures=$((failures + 1))
fi

( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

echo
if (( failures == 0 ))
then
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
