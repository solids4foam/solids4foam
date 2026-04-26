#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# cylindricalPressureVessel regression test
# Checks the final point-displacement magnitude at the inner
# radius probe used by the tutorial.
# ============================================================

DISP_MIN=3.10
DISP_MAX=3.25

ALLRUN_LOGFILE="log.Allrun"

CASES=(
    "displacement::3.10:3.25"
    "pressureDisplacement:pressureDisplacement:2.20:2.32"
    "pressureDisplacementLinear:pressureDisplacementLinear:0.15:0.17"
    "pressureDisplacementUnsteady:pressureDisplacementUnsteady:1.50:1.60"
)

echo "============================================================"
echo "cylindricalPressureVessel regression test"
echo "Final probe displacement magnitude in [${DISP_MIN}, ${DISP_MAX}]"
echo "============================================================"
echo

prepare_case() {
    local case_dir="$1"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done
}

run_case() {
    local case_name="$1"
    local allrun_arg="$2"
    local case_dir="${REGRESSION_ROOT}/${case_name}"

    prepare_case "${case_dir}"
    ( cd "${case_dir}" && ./Allclean > /dev/null 2>&1 ) || true

    if [[ -n "${allrun_arg}" ]]; then
        ( cd "${case_dir}" && ./Allrun "${allrun_arg}" > "${ALLRUN_LOGFILE}" 2>&1 )
    else
        ( cd "${case_dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
    fi

    echo "${case_dir}"
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
    rm -rf "${REGRESSION_ROOT}"
    mkdir -p "${REGRESSION_ROOT}"

    for case_spec in "${CASES[@]}"; do
        IFS=':' read -r case_name allrun_arg min_value max_value \
            <<< "${case_spec}"
        run_case "${case_name}" "${allrun_arg}" > /dev/null
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

failures=0

check_range() {
    local case_name="$1"
    local label="$2"
    local value="$3"
    local min_value="$4"
    local max_value="$5"

    if [[ -z "${value}" ]]; then
        echo "FAIL: ${case_name}: could not extract ${label}"
        failures=$((failures + 1))
    elif awk "BEGIN {exit !(${value} >= ${min_value} && ${value} <= ${max_value})}"; then
        printf "PASS: %s: %s = %.6g\n" "${case_name}" "${label}" "${value}"
    else
        printf "FAIL: %s: %s = %.6g\n" "${case_name}" "${label}" "${value}"
        failures=$((failures + 1))
    fi
}

extract_final_probe_displacement() {
    local case_dir="$1"
    local value_file

    value_file=$(find "${case_dir}/postProcessing" \
        -name 'solidPointDisplacement_*.dat' -print 2>/dev/null \
        | tail -n 1)

    if [[ -z "${value_file}" ]]; then
        return
    fi

    awk 'END {print $5}' "${value_file}"
}

for case_spec in "${CASES[@]}"; do
    IFS=':' read -r case_name allrun_arg min_value max_value \
        <<< "${case_spec}"
    case_dir="${REGRESSION_ROOT}/${case_name}"

    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: ${case_name}"
        continue
    fi

    check_range \
        "${case_name}" "final probe displacement magnitude" \
        "$(extract_final_probe_displacement "${case_dir}")" \
        "${min_value}" "${max_value}"
done

if [ "$CHECK_ONLY" = false ]; then
    for case_dir in "${REGRESSION_ROOT}"/*; do
        if [[ -d "${case_dir}" ]]; then
            ( cd "${case_dir}" && ./Allclean > /dev/null 2>&1 ) || true
        fi
    done
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
