#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# Elastoplastic perforated plate regression test
# Checks strain, stress, and plastic yielding
# ============================================================

# Reference ranges (order-of-magnitude + robustness)
EPSILON_MIN=4.5e-3
EPSILON_MAX=6.0e-3

SIGMA_MIN=1.5e8
SIGMA_MAX=2.2e8

YIELD_MIN=28
YIELD_MAX=32

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Elastoplastic perforated plate regression test"
echo "Max epsilonEq           in [${EPSILON_MIN}, ${EPSILON_MAX}]"
echo "Max sigmaEq (von Mises) in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "Yielding cells          in [${YIELD_MIN}, ${YIELD_MAX}]"
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

prepare_case

# Clean case
( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

# Run case
( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1
}

extract_yielding_cells() {
    grep "cells .* are actively yielding" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 101 \
        | head -n 1 \
        | awk '{print $1}'
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

epsilon=$(extract_max_epsilon)
sigma=$(extract_max_sigma)
yielding_cells=$(extract_yielding_cells)

if [[ -z "${epsilon}" || -z "${sigma}" || -z "${yielding_cells}" ]]
then
    echo "FAIL: Could not extract one or more regression quantities"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

# --- epsilonEq ---
if awk "BEGIN {exit !(${epsilon} >= ${EPSILON_MIN} && ${epsilon} <= ${EPSILON_MAX})}"
then
    printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
else
    printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
    failures=$((failures + 1))
fi

# --- sigmaEq ---
if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${sigma} <= ${SIGMA_MAX})}"
then
    printf "PASS: Max sigmaEq = %.6g\n" "${sigma}"
else
    printf "FAIL: Max sigmaEq = %.6g\n" "${sigma}"
    failures=$((failures + 1))
fi

# --- yielding cells ---
if (( yielding_cells >= YIELD_MIN && yielding_cells <= YIELD_MAX ))
then
    printf "PASS: Yielding cells = %d\n" "${yielding_cells}"
else
    printf "FAIL: Yielding cells = %d\n" "${yielding_cells}"
    failures=$((failures + 1))
fi

# Clean case again
( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

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
