#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Plate-with-hole regression test
# Checks numerical vs analytical solution
# ============================================================

# Regression tolerances
DISP_TOL=1e-7
POINT_DISP_TOL=1e-7
STRESS_TOL=2e5   # component-0 LInf

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Plate-with-hole regression test"
echo "DDifference LInf         < ${DISP_TOL}"
echo "pointDDifference LInf   < ${POINT_DISP_TOL}"
echo "Stress component-0 LInf < ${STRESS_TOL}"
echo "============================================================"
echo

# Clean case
./Allclean > /dev/null 2>&1 || true

# Run case
./Allrun > "${ALLRUN_LOGFILE}" 2>&1

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_disp_linf() {
    local field="$1"
    grep -A2 "Writing ${field} field" "${SOLVER_LOGFILE}" \
        | grep "Norms:" -A1 \
        | tail -n 1 \
        | awk '{print $3}'
}

extract_stress_linf_comp0() {
    grep -A6 "Writing cellStressDifference field" "${SOLVER_LOGFILE}" \
        | awk '
            /Component:[[:space:]]*0/ {getline; getline; print $3}
        '
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

disp_linf=$(extract_disp_linf "DDifference")
point_disp_linf=$(extract_disp_linf "pointDDifference")
stress_linf=$(extract_stress_linf_comp0)

if [[ -z "${disp_linf}" || -z "${point_disp_linf}" || -z "${stress_linf}" ]]; then
    echo "FAIL: Could not extract one or more error norms"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${disp_linf} < ${DISP_TOL})}"; then
    printf "PASS: DDifference LInf = %.6g\n" "${disp_linf}"
else
    printf "FAIL: DDifference LInf = %.6g\n" "${disp_linf}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${point_disp_linf} < ${POINT_DISP_TOL})}"; then
    printf "PASS: pointDDifference LInf = %.6g\n" "${point_disp_linf}"
else
    printf "FAIL: pointDDifference LInf = %.6g\n" "${point_disp_linf}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${stress_linf} < ${STRESS_TOL})}"; then
    printf "PASS: stress component-0 LInf = %.6g\n" "${stress_linf}"
else
    printf "FAIL: stress component-0 LInf = %.6g\n" "${stress_linf}"
    failures=$((failures + 1))
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
