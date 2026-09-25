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
# squarePlate regression test
# Checks the peak transverse deflection written to wVf.
# ============================================================

WF_MIN=6.8e-4
WF_MAX=7.0e-4

# The answer of the removed legacy mechanicalModel, from the last commit that
# had it (mcl-stage8-coverage, c3a92b3d), on OpenFOAM.com v2512, the one fork
# this case runs on. kirchhoffPlate reads rho, E and nu from its single
# linearElastic material rather than asking for a stress, and the deflection is
# linear in those constants, so the framework reproduced this in every written
# digit. The value is written to six figures, so a round-off difference on
# another machine can move its last digit; the tolerance, 1e-5 relative, is
# a few units in that digit, and a 0.001% change in E or nu still reaches it
LEGACY_MAX_WVF=0.000691225
LEGACY_REL_TOL=1e-5

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "squarePlate regression test"
echo "Max wVf in [${WF_MIN}, ${WF_MAX}]"
echo "============================================================"
echo

prepare_case() {
    local d="$1"
    rm -rf "${d}"
    mkdir -p "${d}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
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

field_file=$(find "${CASE_DIR}" -path '*/1/wVf' -print | tail -n 1)
if [[ -z "${field_file}" ]]; then
    echo "FAIL: Could not find wVf field output"
    exit 1
fi

max_wvf=$(awk '
    BEGIN {inlist=0; max=""}
    /^\($/ {inlist=1; next}
    /^\)$/ {inlist=0; next}
    inlist && $1 ~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/ {
        if (max == "" || $1 > max) {
            max = $1
        }
    }
    END {print max}
' "${field_file}")

if [[ -z "${max_wvf}" ]]; then
    echo "FAIL: Could not extract max wVf"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${max_wvf} >= ${WF_MIN} && ${max_wvf} <= ${WF_MAX})}"; then
    printf "PASS: Max wVf = %.6g\n" "${max_wvf}"
else
    printf "FAIL: Max wVf = %.6g\n" "${max_wvf}"
    failures=$((failures + 1))
fi

# ------------------------------------------------------------
# Against the legacy answer
#
# kirchhoffPlate takes rho, E and nu from the mechanicalConstitutiveLaw
# framework, so this checks that it reads the same constants there as it did
# from the legacy model: the same constants give the same number, not merely
# one inside the band
# ------------------------------------------------------------

solver_log=$(find "${CASE_DIR}" -name 'log.solids4Foam' | tail -n 1)

if [[ -n "${solver_log}" ]] \
    && grep -q "Selecting mechanical constitutive law" "${solver_log}"
then
    echo "PASS: the material came from the framework"
else
    echo "FAIL: the solver log shows no mechanical constitutive law"
    failures=$((failures + 1))
fi

if awk "BEGIN {d = ${max_wvf} - ${LEGACY_MAX_WVF}; if (d < 0) d = -d;
               exit !(d <= ${LEGACY_REL_TOL} * ${LEGACY_MAX_WVF})}"
then
    printf "PASS: Max wVf matches the legacy model (%.8g vs %.8g)\n" \
        "${max_wvf}" "${LEGACY_MAX_WVF}"
else
    printf "FAIL: Max wVf differs from the legacy model (%.8g vs %.8g)\n" \
        "${max_wvf}" "${LEGACY_MAX_WVF}"
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
