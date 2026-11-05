#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# healthyHuman regression test
# The left ventricle over the last cycle, [19 T, 20 T] with T = 1 s: the
# end-diastolic and end-systolic volumes and the ejection fraction
# ============================================================

# ------------------------------------------------------------
# Regression tolerances and reference values
# ------------------------------------------------------------

REF_EDV=141.707     # mL
REF_ESV=55.843      # mL
REF_EF=0.6059
VOLUME_TOL=0.05     # mL
EF_TOL=0.001

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
RESULTS_FILE="postProcessing/circulation/0/circulation.dat"

echo "============================================================"
echo "healthyHuman regression test"
echo "EDV, ESV within ${VOLUME_TOL} mL, EF within ${EF_TOL}"
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
            regressionTests|postProcessing|plots|log.*)
                continue
                ;;
        esac
        if [[ "${base_item}" =~ ^[0-9.eE+-]+$ ]]; then
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
# Extract values
# ------------------------------------------------------------

values=""
if [[ -f "${CASE_DIR}/${RESULTS_FILE}" ]]; then
    values=$(awk '
        /^# t/ { for (i = 2; i <= NF; i++) if ($i == "LV.V") c = i - 1 }
        /^#/ || NF == 0 { next }
        c > 0 && $1 >= 19 - 1e-9 && $1 <= 20 + 1e-9 {
            if (n == 0 || $c > hi) hi = $c
            if (n == 0 || $c < lo) lo = $c
            n++
        }
        END { if (n > 0) printf "%.6f %.6f %.6f", hi, lo, (hi - lo)/hi }
    ' "${CASE_DIR}/${RESULTS_FILE}")
fi

if [[ -z "${values}" ]]; then
    echo "FAIL: Could not extract the last cycle of LV.V from ${RESULTS_FILE}"
    echo
    echo "============================================================"
    echo "Regression test FAILED"
    echo "============================================================"
    exit 1
fi

IFS=' ' read -r edv esv ef <<< "${values}"

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check() {
    local name="$1" value="$2" ref="$3" tol="$4"
    if awk -v v="${value}" -v r="${ref}" -v t="${tol}" \
        'BEGIN { d = v - r; exit !((d < 0 ? -d : d) < t) }'
    then
        printf "PASS: %s = %s (reference %s)\n" "${name}" "${value}" "${ref}"
    else
        printf "FAIL: %s = %s (reference %s)\n" "${name}" "${value}" "${ref}"
        failures=$((failures + 1))
    fi
}

check "EDV [mL]" "${edv}" "${REF_EDV}" "${VOLUME_TOL}"
check "ESV [mL]" "${esv}" "${REF_ESV}" "${VOLUME_TOL}"
check "EF" "${ef}" "${REF_EF}" "${EF_TOL}"

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
