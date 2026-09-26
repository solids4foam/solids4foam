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

# The lower bound also guards the hydrostatic stress smoothing that the case
# asks for with solvePressureEqn: with it the probe reads 3.18206, as it did on
# the removed legacy mechanicalModel, and without it 3.16927. The log line
# that says the smoothing is on is checked too, below
DISP_MIN=3.175
DISP_MAX=3.25

ALLRUN_LOGFILE="log.Allrun"

CASES=(
    "displacement::${DISP_MIN}:${DISP_MAX}"
    "pressureDisplacement:pressureDisplacement:2.20:2.32"
    "pressureDisplacementLinear:pressureDisplacementLinear:0.15:0.17"
    "pressureDisplacementUnsteady:pressureDisplacementUnsteady:1.50:1.60"
)

# The final probe displacement of the removed legacy mechanicalModel, from the
# last commit that had it (mcl-stage8-coverage, c3a92b3d), on foam-extend 4.1,
# where coupledPressureDisplacementSolid runs. Each is a case above, as
# name : legacy value : relative tolerance, the tolerance being the one the
# framework was held to against it
LEGACY_REFERENCES=(
    # Nonlinear: the framework takes the law's isochoric stress minus the
    # solved pressure, where the legacy law's pressureDisplacement mode used
    # mu*(b - I)/J, which is not deviatoric. At nu = 0.5 the two differ by
    # about 1e-4 in this probe, so the two are not the same answer and the
    # bound is the 1e-3 that allowed for it
    "pressureDisplacement:2.26015672427:1e-3"
    # Linear: the law is not evaluated and the stiffness is the same shear
    # modulus, so the framework reproduced the legacy value to round-off. This
    # is the legacy value CI logged (foam-extend-4.1-PETSc image,
    # mcl-stage8-coverage 652cdb62); a macOS foam-extend 4.1 build gave
    # 0.159759049445, 3e-3 apart, so it holds for the CI build. The tolerance,
    # 1e-6, is well above round-off and below that difference
    "pressureDisplacementLinear:0.160195657336:1e-6"
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

# The displacement arm asks for the hydrostatic stress smoothing, and the
# solid model has to say it is doing it: the band above is set to catch its
# loss, and this is the direct evidence
displacement_dir="${REGRESSION_ROOT}/displacement"
if ! solids4Foam::regressionCaseSkipped "${displacement_dir}/${ALLRUN_LOGFILE}"
then
    if grep -q "smoothing the hydrostatic stress (solvePressureEqn)" \
        "${displacement_dir}/log.solids4Foam" 2>/dev/null
    then
        echo "PASS: displacement: the hydrostatic stress is smoothed"
    else
        echo "FAIL: displacement: no hydrostatic stress smoothing in the log"
        failures=$((failures + 1))
    fi
fi

# The pressure-displacement arms against the legacy answers
for reference in "${LEGACY_REFERENCES[@]}"; do
    IFS=':' read -r case_name legacy_value rel_tol <<< "${reference}"
    case_dir="${REGRESSION_ROOT}/${case_name}"

    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: ${case_name} against the legacy model"
        continue
    fi

    if ! grep -q "taking the stiffness from the mechanicalConstitutiveLaw" \
        "${case_dir}/log.solids4Foam" 2>/dev/null
    then
        echo "FAIL: ${case_name}: did not take its stiffness from the framework"
        failures=$((failures + 1))
        continue
    fi

    value="$(extract_final_probe_displacement "${case_dir}")"

    if [[ -z "${value}" ]]; then
        echo "FAIL: ${case_name}: could not extract the probe value"
        failures=$((failures + 1))
    elif awk "BEGIN {d = ${value} - ${legacy_value}; \
        exit !(${value} > 0 && d*d <= (${rel_tol}*${legacy_value})^2)}"
    then
        printf "PASS: %s against the legacy model: %.12g vs %.12g\n" \
            "${case_name}" "${value}" "${legacy_value}"
    else
        printf "FAIL: %s against the legacy model: %.12g vs %.12g (tolerance %s)\n" \
            "${case_name}" "${value}" "${legacy_value}" "${rel_tol}"
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
