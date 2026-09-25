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
    "pressureDisplacementManager:pressureDisplacement:2.20:2.32"
    "pressureDisplacementLinearManager:pressureDisplacementLinear:0.15:0.17"
)

# The *Manager arms run the same case with coupledPressureDisplacementSolid
# taking its constitutive response from the mechanicalConstitutiveLaw
# framework, and are compared with their legacy twin below. The switch goes in
# the solid model's coeffs sub-dictionary, which is where the model reads it
FRAMEWORK_PAIRS=(
    # legacy arm : framework arm : relative tolerance
    #
    # Nonlinear: the framework takes the law's isochoric stress minus the
    # solved pressure, where the legacy law's pressureDisplacement mode uses
    # mu*(b - I)/J, which is not deviatoric. At nu = 0.5 the two differ by
    # about 1e-4 in this probe
    "pressureDisplacement:pressureDisplacementManager:1e-3"
    # Linear: the law is not evaluated and the stiffness is the same shear
    # modulus, so the two agree to round-off
    "pressureDisplacementLinear:pressureDisplacementLinearManager:1e-9"
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

    if [[ "${case_name}" == *Manager ]]; then
        local dict
        local switch="    useMechanicalConstitutiveLawManager yes;"
        for dict in \
            "${case_dir}"/caseOptions/pressureDisplacement*/*/constant/solidProperties
        do
            sed -i \
                "/coupledPressureDisplacementSolidCoeffs/,/{/ s|{|{\n${switch}|" \
                "${dict}"
        done
    fi

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

# Each framework arm against its legacy twin. Asserted positively that each
# took the path it was set up for, or the comparison is a run against itself
for pair in "${FRAMEWORK_PAIRS[@]}"; do
    IFS=':' read -r legacy_name framework_name rel_tol <<< "${pair}"
    legacy_dir="${REGRESSION_ROOT}/${legacy_name}"
    framework_dir="${REGRESSION_ROOT}/${framework_name}"

    if solids4Foam::regressionCaseSkipped "${framework_dir}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: ${framework_name} against ${legacy_name}"
        continue
    fi

    if ! grep -q "taking the stiffness from the mechanicalConstitutiveLaw" \
        "${framework_dir}/log.solids4Foam" 2>/dev/null
    then
        echo "FAIL: ${framework_name}: did not use the framework"
        failures=$((failures + 1))
        continue
    fi

    if grep -q "Creating the mechanicalConstitutiveLawManager" \
        "${legacy_dir}/log.solids4Foam" 2>/dev/null
    then
        echo "FAIL: ${legacy_name}: used the framework"
        failures=$((failures + 1))
        continue
    fi

    legacy_value="$(extract_final_probe_displacement "${legacy_dir}")"
    framework_value="$(extract_final_probe_displacement "${framework_dir}")"

    if [[ -z "${legacy_value}" || -z "${framework_value}" ]]; then
        echo "FAIL: ${framework_name}: could not extract both probe values"
        failures=$((failures + 1))
    elif awk "BEGIN {d = ${framework_value} - ${legacy_value}; \
        exit !(${legacy_value} > 0 && d*d <= (${rel_tol}*${legacy_value})^2)}"
    then
        printf "PASS: %s against %s: %.12g vs %.12g\n" \
            "${framework_name}" "${legacy_name}" \
            "${framework_value}" "${legacy_value}"
    else
        printf "FAIL: %s against %s: %.12g vs %.12g (tolerance %s)\n" \
            "${framework_name}" "${legacy_name}" \
            "${framework_value}" "${legacy_value}" "${rel_tol}"
        failures=$((failures + 1))
    fi
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
