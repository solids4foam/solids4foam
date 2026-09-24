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
# curvedBeams regression test
# Uses the reaction force history as a cheap contact benchmark check.
# ============================================================

FORCE_Y_MIN=-17.8
FORCE_Y_MAX=-17.6

# The final force_y of the removed legacy mechanicalModel, from the last
# commit that had it (mcl-stage8-coverage, c3a92b3d), on foam-extend 4.1, the
# one fork this case runs on. The band above is a correctness bound; this is a
# much tighter one, because the framework solves the same problem. It
# reproduced this value to the eight digits printed; the tolerance is the
# 0.02 N the two arms were held to when both ran
LEGACY_FINAL_FORCE_Y=-17.724965
LEGACY_FORCE_TOL=0.02

ALLRUN_LOGFILE="log.Allrun"
FORCE_FILE="postProcessing/0/solidForcesdisplacement.dat"

echo "============================================================"
echo "curvedBeams regression test"
echo "Final displacement-patch force_y in [${FORCE_Y_MIN}, ${FORCE_Y_MAX}] N"
echo "============================================================"
echo

prepare_case() {
    local d="${1:-${CASE_DIR}}"
    rm -rf "${d}"
    mkdir -p "${d}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${d}/"
    done

    sed -i.bak 's/^endTime[[:space:]]\+31.5;/endTime         9;/' "${d}/system/controlDict"
    rm -f "${d}/system/controlDict.bak"
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

if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
    echo "FAIL: Could not find ${FORCE_FILE}"
    exit 1
fi

final_force_y=$(awk 'END {print $3}' "${CASE_DIR}/${FORCE_FILE}")

if [[ -z "${final_force_y}" ]]; then
    echo "FAIL: Could not extract final force_y"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${final_force_y} >= ${FORCE_Y_MIN} && ${final_force_y} <= ${FORCE_Y_MAX})}"; then
    printf "PASS: Final force_y = %.6g\n" "${final_force_y}"
else
    printf "FAIL: Final force_y = %.6g\n" "${final_force_y}"
    failures=$((failures + 1))
fi

# ------------------------------------------------------------
# Against the legacy answer
#
# This is the framework's coverage of a contact case. It matters beyond this
# tutorial: the contact penalty models look impK up from the registry by name,
# so frameworkImpK() has to register a field of the same name, with the same
# dimensions and boundary types, as the legacy impK() it replaced. Nothing else
# tests that, and the same lookup is used by the cohesive zone models and by
# elasticWallPressure in FSI.
#
# The material is also history dependent, so this exercises a framework law
# carrying plastic state through a contact solve
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

# The band is wide enough to hold answers that disagree materially, so the
# legacy answer is checked too
if awk "BEGIN {d = ${final_force_y} - ${LEGACY_FINAL_FORCE_Y};
               if (d < 0) d = -d;
               exit !(d <= ${LEGACY_FORCE_TOL})}"
then
    printf "PASS: force_y matches the legacy model (%.8g vs %.8g)\n" \
        "${final_force_y}" "${LEGACY_FINAL_FORCE_Y}"
else
    printf "FAIL: force_y differs from the legacy model (%.8g vs %.8g)\n" \
        "${final_force_y}" "${LEGACY_FINAL_FORCE_Y}"
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
