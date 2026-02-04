#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Contact patch test regression test
# Checks numerical vs analytical solution
# ============================================================

# Reference ranges
EPSILON_MIN=0.0065
EPSILON_MAX=0.0067

SIGMA_MIN=9800
SIGMA_MAX=10000

SIGMA_Y_REL_ERROR_MIN=0
SIGMA_Y_REL_ERROR_MAX=2

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Contact patch test regression test"
echo "Max epsilonEq           in [${EPSILON_MIN}, ${EPSILON_MAX}]"
echo "Max sigmaEq (von Mises) in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "Max sigma_y relative error (in %) in " \
     "[${SIGMA_Y_REL_ERROR_MIN}, ${SIGMA_Y_REL_ERROR_MAX}]"
echo "============================================================"
echo

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
    ./Allclean > /dev/null 2>&1 || true
    ./Allrun > "${ALLRUN_LOGFILE}" 2>&1
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_epsilon() {
    grep "Max epsilonEq" "${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1
}

extract_relative_sigma_y_error() {
    grep "Average relative error in sigma_y field:" "${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tr -d '%' \
        | tail -n 1
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

epsilon=$(extract_max_epsilon)
sigma=$(extract_max_sigma)
sigma_y_rel_error=$(extract_relative_sigma_y_error)

if [[ -z "${epsilon}" || -z "${sigma}" || -z "${sigma_y_rel_error}" ]]
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

# --- average relative sigma_y error (percent) ---
if awk "BEGIN {exit !(${sigma_y_rel_error} >= ${SIGMA_Y_REL_ERROR_MIN} && \
${sigma_y_rel_error} <= ${SIGMA_Y_REL_ERROR_MAX})}"
then
    printf "PASS: Avg rel. error sigma_y = %.6g%%\n" \
        "${sigma_y_rel_error}"
else
    printf "FAIL: Avg rel. error sigma_y = %.6g%%\n" \
        "${sigma_y_rel_error}"
    failures=$((failures + 1))
fi

echo
if (( failures == 0 ));
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
