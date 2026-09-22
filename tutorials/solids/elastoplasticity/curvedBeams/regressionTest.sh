#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
FRAMEWORK_DIR="${REGRESSION_ROOT}/framework"
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

# Read the final displacement-patch reaction force from a case
read_final_force_y() {
    local d="$1"
    [[ -f "${d}/${FORCE_FILE}" ]] || return 1
    awk 'END {print $3}' "${d}/${FORCE_FILE}"
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

    # The framework arm differs in this one entry and nothing else. It is
    # applied after Allclean, which restores the stored dictionaries
    prepare_case "${FRAMEWORK_DIR}"
    ( cd "${FRAMEWORK_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    sed -i.bak \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${FRAMEWORK_DIR}/constant/solidProperties"
    rm -f "${FRAMEWORK_DIR}/constant/solidProperties.bak"

    if ! grep -q "useMechanicalConstitutiveLawManager" \
        "${FRAMEWORK_DIR}/constant/solidProperties"
    then
        echo "FAIL: could not set the framework switch on the framework arm"
        exit 1
    fi

    ( cd "${FRAMEWORK_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
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
# The framework arm
#
# This is the first framework coverage of a contact case. It matters beyond
# this tutorial: the contact penalty models look impK up from the registry by
# name, so frameworkImpK() has to register a field of the same name, with the
# same dimensions and boundary types, as the legacy impK() it replaces.
# Nothing else tests that, and the same lookup is used by the cohesive zone
# models and by elasticWallPressure in FSI.
#
# The material is also history dependent, so this exercises a framework law
# carrying plastic state through a contact solve
# ------------------------------------------------------------

if [ "$CHECK_ONLY" = false ]; then
    if solids4Foam::regressionCaseSkipped "${FRAMEWORK_DIR}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: framework arm skipped in this environment"
    else
        fw_log=$(find "${FRAMEWORK_DIR}" -name 'log.solids4Foam' | tail -n 1)
        main_log=$(find "${CASE_DIR}" -name 'log.solids4Foam' | tail -n 1)

        if [[ -n "${main_log}" ]] \
            && grep -q "mechanicalConstitutiveLawManager" "${main_log}"
        then
            echo "FAIL: the legacy arm used the framework"
            failures=$((failures + 1))
        else
            echo "PASS: legacy arm took the legacy path"
        fi

        if [[ -n "${fw_log}" ]] \
            && grep -q "mechanicalConstitutiveLawManager" "${fw_log}"
        then
            echo "PASS: framework arm took the framework path"
        else
            echo "FAIL: framework arm did not take the framework path"
            failures=$((failures + 1))
        fi

        fw_force_y=$(read_final_force_y "${FRAMEWORK_DIR}" || true)

        if [[ -z "${fw_force_y}" ]]; then
            echo "FAIL: framework arm produced no force history"
            failures=$((failures + 1))
        elif awk "BEGIN {exit !(${fw_force_y} >= ${FORCE_Y_MIN} \
                  && ${fw_force_y} <= ${FORCE_Y_MAX})}"
        then
            printf "PASS: framework final force_y = %.6g\n" "${fw_force_y}"
        else
            printf "FAIL: framework final force_y = %.6g\n" "${fw_force_y}"
            failures=$((failures + 1))
        fi
    fi
fi

if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${FRAMEWORK_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
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
