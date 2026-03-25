#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

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

# Clean case again
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
