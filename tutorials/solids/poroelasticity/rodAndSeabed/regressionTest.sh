#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# Rod and seabed regression test
# Checks strain and stress
# ============================================================

# Reference ranges (order-of-magnitude + robustness)
EPSILON_MIN=6e-5
EPSILON_MAX=1.1e-4

SIGMA_MIN=70e3
SIGMA_MAX=90e3

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Road and seabed regression test"
echo "Max epsilonEq           in [${EPSILON_MIN}, ${EPSILON_MAX}]"
echo "Max sigmaEq (von Mises) in [${SIGMA_MIN}, ${SIGMA_MAX}]"
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

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

epsilon=$(extract_max_epsilon)
sigma=$(extract_max_sigma)

if [[ -z "${epsilon}" || -z "${sigma}" ]]
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
