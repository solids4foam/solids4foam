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
# annularPlate regression test
# Checks the transverse deflection on the free outer edge of the default
# clampedInner configuration, written to the finite area field w.
# ============================================================

# The 16 x 152 mesh overpredicts the exact free-edge deflection of 7.503e-6 m
# by about 0.7%; the band allows a 1% drift either side of the recorded value
WF_MIN=7.48e-6
WF_MAX=7.63e-6

# The value recorded on OpenFOAM.com v2606, written to eight figures. The
# segregated solver converges to a relative tolerance of 1e-8, so a round-off
# difference on another machine moves the last figures only
REF_MAX_W=7.5544528e-06
REF_REL_TOL=1e-4

PATCH="outer"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "annularPlate regression test"
echo "Max w on patch ${PATCH} in [${WF_MIN}, ${WF_MAX}]"
echo "============================================================"
echo

prepare_case() {
    local d="$1"
    rm -rf "${d}"
    mkdir -p "${d}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" || "${base_item}" == "verification" ]]; then
            continue
        fi
        cp -a "${item}" "${d}/"
    done
}

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
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

field_file="${CASE_DIR}/1/finite-area/w"
if [[ ! -f "${field_file}" ]]; then
    echo "FAIL: Could not find the finite area field w"
    exit 1
fi

# The largest value in the nonuniform value list of the patch
max_w=$(awk -v patch="${PATCH}" '
    BEGIN {inpatch=0; inlist=0; max=""}
    $1 == patch {inpatch=1; next}
    inpatch && /^\s*}/ {inpatch=0; inlist=0; next}
    inpatch && /^\($/ {inlist=1; next}
    inpatch && /^\)$/ {inlist=0; next}
    inpatch && inlist && $1 ~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/ {
        if (max == "" || $1 > max) {
            max = $1
        }
    }
    END {print max}
' "${field_file}")

if [[ -z "${max_w}" ]]; then
    echo "FAIL: Could not extract max w on patch ${PATCH}"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${max_w} >= ${WF_MIN} && ${max_w} <= ${WF_MAX})}"; then
    printf "PASS: Max w on patch %s = %.6g\n" "${PATCH}" "${max_w}"
else
    printf "FAIL: Max w on patch %s = %.6g\n" "${PATCH}" "${max_w}"
    failures=$((failures + 1))
fi

if awk "BEGIN {d = ${max_w} - ${REF_MAX_W}; if (d < 0) d = -d;
               exit !(d <= ${REF_REL_TOL} * ${REF_MAX_W})}"
then
    printf "PASS: Max w matches the reference (%.8g vs %.8g)\n" \
        "${max_w}" "${REF_MAX_W}"
else
    printf "FAIL: Max w differs from the reference (%.8g vs %.8g)\n" \
        "${max_w}" "${REF_MAX_W}"
    failures=$((failures + 1))
fi

solver_log=$(find "${CASE_DIR}" -name 'log.solids4Foam' | tail -n 1)

if [[ -n "${solver_log}" ]] \
    && grep -q "The residuals have converged" "${solver_log}"
then
    echo "PASS: the correctors converged"
else
    echo "FAIL: the solver log shows no converged correctors"
    failures=$((failures + 1))
fi

if [ "$CHECK_ONLY" = false ]; then
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
