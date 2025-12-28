#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# viscoTube regression test
# Checks order of magnitude of strain and von Mises stress
# ============================================================

# Order-of-magnitude tolerances
EPS_MIN=1e-5
EPS_MAX=1e-3

SIGMA_MIN=1e6
SIGMA_MAX=1e8

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "viscoTube regression test"
echo "epsilonEq expected: ${EPS_MIN} < eps < ${EPS_MAX}"
echo "sigmaEq expected  : ${SIGMA_MIN} < sigma < ${SIGMA_MAX}"
echo "============================================================"
echo

# Clean case
./Allclean > /dev/null 2>&1 || true

# Run case
./Allrun > "${ALLRUN_LOGFILE}" 2>&1

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_epsilon() {
    grep "Max epsilonEq" "${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk -F '=' '{print $2}' \
        | tr -d '[:space:]'
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${SOLVER_LOGFILE}" \
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

if awk "BEGIN {exit !(${epsilon} > ${EPS_MIN} && ${epsilon} < ${EPS_MAX})}"
then
    printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
else
    printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${sigma} > ${SIGMA_MIN} && ${sigma} < ${SIGMA_MAX})}"
then
    printf "PASS: Max sigmaEq = %.6g\n" "${sigma}"
else
    printf "FAIL: Max sigmaEq = %.6g\n" "${sigma}"
    failures=$((failures + 1))
fi

# Clean case again
./Allclean > /dev/null 2>&1 || true

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
