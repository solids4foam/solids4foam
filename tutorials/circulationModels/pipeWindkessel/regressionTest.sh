#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
RESTART_DIR="${REGRESSION_ROOT}/restart"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# pipeWindkessel regression test
# 1. The outlet pressure against the exact solution of the Windkessel
# 2. A restart: the case run to t = 0.5 on two ranks, reconstructed, and
#    continued in serial to t = 1, against the continuous run of 1.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances and reference values
# ------------------------------------------------------------

# RMSE of the outlet pressure Po over 0 < t <= 1 against the exact solution
# (see analyse.py); about 1 % of the amplitude of Po
REF_RMSE=0.3385
RMSE_TOL=0.005

# Largest difference of a restarted row from the continuous one, relative to
# the range of its column: round-off, the restart being exact
RESTART_TOL=1e-9

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
COUPLING_FILE="postProcessing/circulationCoupling/0/coupling.dat"
RESTART_FILE="postProcessing/circulationCoupling/0.5/coupling.dat"

echo "============================================================"
echo "pipeWindkessel regression test"
echo "|RMSE - ${REF_RMSE}| < ${RMSE_TOL}, restart difference < ${RESTART_TOL}"
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

set_entry() {
    foamDictionary -entry "$2" -set "$3" "$1" > /dev/null
}

# ------------------------------------------------------------
# Clean & run cases
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

if [ "$CHECK_ONLY" = false ]; then
    prepare_case "${RESTART_DIR}"
    (
        cd "${RESTART_DIR}"
        set_entry system/decomposeParDict numberOfSubdomains 2
        set_entry system/decomposeParDict method scotch
        set_entry system/controlDict endTime 0.5
        set_entry system/controlDict writeInterval 10
        solids4Foam::runApplication blockMesh
        solids4Foam::runApplication decomposePar
        solids4Foam::runParallel pimpleFoam
        solids4Foam::runApplication reconstructPar -latestTime
        rm -rf processor*
        set_entry system/controlDict endTime 1
        set_entry system/controlDict startFrom latestTime
        solids4Foam::runApplication -s restart pimpleFoam
    ) > "${RESTART_DIR}/${ALLRUN_LOGFILE}" 2>&1 || true
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

# 1. The outlet pressure against the exact solution
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

# 2. The restart: the rows of (0.5, 1] against the continuous run
if grep -q "circulationCoupler: resumed from" \
    "${RESTART_DIR}/log.pimpleFoam.restart" 2> /dev/null
then
    echo "PASS: the restarted run resumed from the circulation checkpoint"
else
    echo "FAIL: the restarted run did not resume from the circulation checkpoint"
    failures=$((failures + 1))
fi

if [[ -f "${RESTART_DIR}/${RESTART_FILE}" ]]; then
    restart_diff=$(awk '
        /^#/ || NF == 0 { next }
        FNR == NR {
            if ($1 > 0.5 + 1e-9) {
                for (i = 1; i <= NF; i++) {
                    ref[$1, i] = $i
                    if (!(i in lo) || $i < lo[i]) lo[i] = $i
                    if (!(i in hi) || $i > hi[i]) hi[i] = $i
                }
                nref++
            }
            next
        }
        {
            if (!(($1, 1) in ref)) { missing++; next }
            for (i = 2; i <= NF; i++) {
                range = hi[i] - lo[i]
                d = $i - ref[$1, i]
                d = (d < 0 ? -d : d)/(range > 0 ? range : 1)
                if (d > worst) worst = d
            }
            n++
        }
        END {
            if (missing > 0 || n == 0 || n != nref) print "rows"
            else printf "%.3g", worst + 0
        }
    ' "${CASE_DIR}/${COUPLING_FILE}" "${RESTART_DIR}/${RESTART_FILE}")
else
    restart_diff="none"
fi

if [[ "${restart_diff}" == "none" || "${restart_diff}" == "rows" ]]; then
    echo "FAIL: the restarted run's rows do not match those of the continuous run"
    failures=$((failures + 1))
elif awk -v d="${restart_diff}" -v tol="${RESTART_TOL}" \
    'BEGIN { exit !(d < tol) }'
then
    printf "PASS: restarted rows differ from the continuous run by %s of range\n" \
        "${restart_diff}"
else
    printf "FAIL: restarted rows differ from the continuous run by %s of range\n" \
        "${restart_diff}"
    failures=$((failures + 1))
fi

# Clean cases again, keeping them when a check failed
if [[ "$CHECK_ONLY" = false && ${failures} -eq 0 ]]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${RESTART_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
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
