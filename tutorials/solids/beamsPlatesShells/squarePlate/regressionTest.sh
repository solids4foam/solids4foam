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
# squarePlate regression test
# Checks the peak transverse deflection written to wVf.
# ============================================================

WF_MIN=6.8e-4
WF_MAX=7.0e-4

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

# Read the peak transverse deflection from a case's wVf field
read_max_wvf() {
    local d="$1"
    local f
    f=$(find "${d}" -path '*/1/wVf' -print | tail -n 1)
    if [[ -z "${f}" ]]; then
        return 1
    fi

    awk '
        BEGIN {inlist=0; max=""}
        /^\($/ {inlist=1; next}
        /^\)$/ {inlist=0; next}
        inlist && $1 ~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/ {
            if (max == "" || $1 > max) {
                max = $1
            }
        }
        END {print max}
    ' "${f}"
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

    # The framework arm differs in this one entry and nothing else. It goes
    # inside the model's coeffs block, which is where the solid model looks
    # for it; at the top level it is read by nothing and silently ignored.
    #
    # After Allclean, not before: Allclean ends in restoreCaseFormat, which
    # puts the stored dictionaries back and would undo this
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
# The framework arm
#
# kirchhoffPlate reads rho, E and nu from its single linearElastic material
# rather than asking for a stress, so this checks that it reads the same
# constants from the constitutive law framework as from the legacy model. The
# plate deflection is linear in those constants, so an equal wVf is an equal
# bending stiffness
# ------------------------------------------------------------

if solids4Foam::regressionCaseSkipped "${FRAMEWORK_DIR}/${ALLRUN_LOGFILE}"; then
    echo "SKIP: framework arm skipped in this environment"
else
    fw_solver_log=$(find "${FRAMEWORK_DIR}" -name 'log.solids4Foam' | tail -n 1)
    main_solver_log=$(find "${CASE_DIR}" -name 'log.solids4Foam' | tail -n 1)

    # Each arm must have taken the path it was set up for, or the comparison
    # is between two copies of the same thing and proves nothing
    if [[ -n "${main_solver_log}" ]] \
        && grep -q "mechanicalConstitutiveLawManager" "${main_solver_log}"
    then
        echo "FAIL: the legacy arm used the framework"
        failures=$((failures + 1))
    else
        echo "PASS: legacy arm took the legacy path"
    fi

    if [[ -n "${fw_solver_log}" ]] \
        && grep -q "mechanicalConstitutiveLawManager" "${fw_solver_log}"
    then
        echo "PASS: framework arm took the framework path"
    else
        echo "FAIL: framework arm did not take the framework path"
        failures=$((failures + 1))
    fi

    max_wvf_fw=$(read_max_wvf "${FRAMEWORK_DIR}" || true)

    if [[ -z "${max_wvf_fw}" ]]; then
        echo "FAIL: framework arm produced no wVf field"
        failures=$((failures + 1))
    elif awk "BEGIN {exit !(${max_wvf_fw} >= ${WF_MIN} && ${max_wvf_fw} <= ${WF_MAX})}"
    then
        printf "PASS: framework Max wVf = %.6g\n" "${max_wvf_fw}"

        # Same constants read two ways, so the answers are the same number,
        # not merely both inside the band
        if awk "BEGIN {d = ${max_wvf_fw} - ${max_wvf}; if (d < 0) d = -d;
                       exit !(d <= 1e-12 * ${max_wvf})}"
        then
            printf "PASS: legacy and framework agree (%.8g vs %.8g)\n" \
                "${max_wvf}" "${max_wvf_fw}"
        else
            printf "FAIL: legacy and framework differ (%.8g vs %.8g)\n" \
                "${max_wvf}" "${max_wvf_fw}"
            failures=$((failures + 1))
        fi
    else
        printf "FAIL: framework Max wVf = %.6g\n" "${max_wvf_fw}"
        failures=$((failures + 1))
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
