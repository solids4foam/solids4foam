#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# pipeWindkesselSolids4Foam regression test
# The outlet pressure against the exact solution of the Windkessel, with the
# pipe solved by solids4Foam and the pimpleFluid fluid model
# ============================================================

# ------------------------------------------------------------
# Regression tolerances and reference values
# ------------------------------------------------------------

# RMSE of the outlet pressure Po over 0 < t <= 1 against the exact solution
# (see ../pipeWindkessel/analyse.py), as for pipeWindkessel
REF_RMSE=0.3385
RMSE_TOL=0.005

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
COUPLING_FILE="postProcessing/circulationCoupling/0/coupling.dat"

echo "============================================================"
echo "pipeWindkesselSolids4Foam regression test"
echo "|RMSE - ${REF_RMSE}| < ${RMSE_TOL}"
echo "============================================================"
echo

prepare_case() {
    local d="$1"
    rm -rf "${d}"
    mkdir -p "${d}"

    # Copy the inputs only: no results of a previous run of the tutorial
    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        case "${base_item}" in
            regressionTests|postProcessing|processor*|log.*)
                continue
                ;;
        esac
        if [[ "${base_item}" != "0" \
           && "${base_item}" =~ ^[0-9.eE+-]+$ ]]; then
            continue
        fi
        cp -a "${item}" "${d}/"
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
    prepare_case "${CASE_DIR}"
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

# The outlet pressure against the exact solution
#    Pc(t) = A R/(1+(wRC)^2) [sin wt - wRC cos wt + wRC exp(-t/RC)]
#    Po(t) = Pc(t) + Rch A sin wt,  A = 1, w = 2 pi, R = 10, Rch = 1, C = 0.01
if [[ -f "${CASE_DIR}/${COUPLING_FILE}" ]]; then
    rmse=$(awk '
        /^# time/ { for (i = 3; i <= NF; i++) if ($i == "outlet.Ppatch") c = i - 1 }
        /^#/ || NF == 0 { next }
        $1 > 0 && $1 <= 1 + 1e-9 {
            A = 1; w = 8*atan2(1, 1); R = 10; Rch = 1; C = 0.01; k = w*R*C
            t = $1
            pc = A*R/(1 + k*k)*(sin(w*t) - k*cos(w*t) + k*exp(-t/(R*C)))
            e = $c - (pc + Rch*A*sin(w*t))
            s += e*e; n++
        }
        END { if (n > 0 && c > 0) printf "%.6g", sqrt(s/n) }
    ' "${CASE_DIR}/${COUPLING_FILE}")
else
    rmse=""
fi

if [[ -z "${rmse}" ]]; then
    echo "FAIL: Could not compute the RMSE from ${COUPLING_FILE}"
    failures=$((failures + 1))
elif awk -v r="${rmse}" -v ref="${REF_RMSE}" -v tol="${RMSE_TOL}" \
    'BEGIN { d = r - ref; exit !((d < 0 ? -d : d) < tol) }'
then
    printf "PASS: RMSE of the outlet pressure = %s (reference %s)\n" \
        "${rmse}" "${REF_RMSE}"
else
    printf "FAIL: RMSE of the outlet pressure = %s (reference %s)\n" \
        "${rmse}" "${REF_RMSE}"
    failures=$((failures + 1))
fi

# Clean case again, keeping it when a check failed
if [[ "$CHECK_ONLY" = false && ${failures} -eq 0 ]]; then
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
