#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"
SOLIDS4FOAM_ROOT_ABS=$(cd "${SCRIPT_DIR}/../../../../" && pwd)

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# Plate-with-hole regression tests
# Checks numerical vs analytical solution for the displacement
# (segregated/segregatedManager/petscSnes/petscSnesPressure/high-order
# variants) and pressure-displacement
#
# Note that petscSnesPressure and pressureDisplacement* are different things:
# petscSnesPressure is the mixed displacement-pressure form of
# linearGeometryTotalDisplacement (solvePressure yes), whereas
# pressureDisplacement* selects coupledPressureDisplacementSolid, which runs
# on foam-extend only.
# solution options.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

DISP_TOL=1e-7
POINT_DISP_TOL=1e-7
STRESS_TOL=2e5
FRAMEWORK_FIELD_ABS_TOL=2e-12

PD_DISP_TOL=3.0e-4
PD_POINT_DISP_TOL=3.0e-4
PD_SIGMA_TOL=2.0e5
PD_P_TOL=9.0e4

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
PARALLEL_N_PROCS=2

APPROACHES=(
    segregated
    segregatedManager
    petscSnes
    petscSnesPressure
    petscSnesPressureManager
    highOrder-movingLeastSquares
    highOrder-kExactLeastSquares
    highOrder-movingLeastSquares-parallel
    highOrder-kExactLeastSquares-parallel
    highOrderFourthOrder
)

PRESSURE_DISPLACEMENT_CASES=(
    "pressureDisplacementCompressible coarse"
    "pressureDisplacementCompressible medium"
    "pressureDisplacementIncompressible coarse"
    "pressureDisplacementIncompressible medium"
    "pressureDisplacementCompressibleManager coarse"
    "pressureDisplacementIncompressibleManager coarse"
)

# Each *Manager arm is its twin run with coupledPressureDisplacementSolid
# taking its constitutive response from the mechanicalConstitutiveLaw
# framework. These cases run the model in linear mode, where the law is not
# evaluated and the stiffness is the same shear modulus on both paths, so the
# error measures must agree with the twin's to round-off
PD_FRAMEWORK_REL_TOL=1e-6

echo "============================================================"
echo "Plate-with-hole regression tests"
echo "DDifference LInf        < ${DISP_TOL}"
echo "pointDDifference LInf   < ${POINT_DISP_TOL}"
echo "Stress component-0 LInf < ${STRESS_TOL}"
echo "Pressure-displacement DError max      < ${PD_DISP_TOL}"
echo "Pressure-displacement pointDError max < ${PD_POINT_DISP_TOL}"
echo "Pressure-displacement sigma*Err max   < ${PD_SIGMA_TOL}"
echo "Pressure-displacement pErr max        < ${PD_P_TOL}"
echo "============================================================"
echo

prepare_case() {
    local case_dir="$1"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        local base_item
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done

    # The regression copy lives deeper than the source tutorial, so the
    # relative SOLIDS4FOAM_ROOT in this local library build no longer points to
    # the repository root.
    if [[ -f "${case_dir}/src/Make/options" ]]; then
        sed -i.bak \
            "s|^SOLIDS4FOAM_ROOT := .*|SOLIDS4FOAM_ROOT := ${SOLIDS4FOAM_ROOT_ABS}|" \
            "${case_dir}/src/Make/options"
    fi

    sed -E -i.bak \
        "s/^[[:space:]]*numberOfSubdomains[[:space:]]+[0-9]+;/numberOfSubdomains ${PARALLEL_N_PROCS};/" \
        "${case_dir}/system/decomposeParDict"
    rm -f "${case_dir}/system/decomposeParDict.bak"

}

run_case() {
    local case_name="$1"
    shift

    local case_dir="${REGRESSION_ROOT}/${case_name}"
    local requested="$1"

    prepare_case "${case_dir}"

    case "${requested}" in
        highOrder-movingLeastSquares|highOrder-kExactLeastSquares)
            local least_squares_type="${requested#highOrder-}"
            sed -E -i.bak \
                "s/^([[:space:]]*)type[[:space:]]+(movingLeastSquares|kExactLeastSquares);/\\1type ${least_squares_type};/" \
                "${case_dir}/constant/solidProperties.highOrder"
            rm -f "${case_dir}/constant/solidProperties.highOrder.bak"
            set -- highOrder "${@:2}"
            ;;
        highOrder-movingLeastSquares-parallel|highOrder-kExactLeastSquares-parallel)
            local least_squares_type="${requested#highOrder-}"
            least_squares_type="${least_squares_type%-parallel}"
            sed -E -i.bak \
                -e "s/^([[:space:]]*)type[[:space:]]+(movingLeastSquares|kExactLeastSquares);/\\1type ${least_squares_type};/" \
                -e "s/^([[:space:]]*)highOrderJacobian[[:space:]]+(true|false);/\\1highOrderJacobian false;/" \
                "${case_dir}/constant/solidProperties.highOrder"
            rm -f "${case_dir}/constant/solidProperties.highOrder.bak"
            set -- highOrder parallel
            ;;
        pressureDisplacement*Manager)
            local switch="    useMechanicalConstitutiveLawManager yes;"
            sed -i \
                "/coupledPressureDisplacementSolidCoeffs/,/{/ s|{|{\n${switch}|" \
                "${case_dir}/caseOptions/pressureDisplacement/hex/common/constant/solidProperties"
            set -- "${requested%Manager}" "${@:2}"
            ;;
    esac

    ( cd "${case_dir}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${case_dir}" && ./Allrun "$@" > "${ALLRUN_LOGFILE}" 2>&1 )

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

    for approach in "${APPROACHES[@]}"; do
        run_case "${approach}" "${approach}" > /dev/null
    done

    for case_args in "${PRESSURE_DISPLACEMENT_CASES[@]}"; do
        IFS=' ' read -r approach mesh <<< "${case_args}"
        run_case "${approach}-${mesh}" "${approach}" "${mesh}" > /dev/null
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

extract_disp_linf() {
    local case_dir="$1"
    local field="$2"

    grep -A2 "Writing ${field} field" "${case_dir}/${SOLVER_LOGFILE}" \
        | grep "Norms:" -A1 \
        | tail -n 1 \
        | awk '{print $3}' \
        || true
}

extract_stress_linf_comp0() {
    local case_dir="$1"

    grep -A6 "Writing cellStressDifference field" \
        "${case_dir}/${SOLVER_LOGFILE}" \
        | tail -6 \
        | awk '
            /Component:[[:space:]]*0/ {getline; getline; print $3}
        ' \
        || true
}

extract_log_value() {
    local case_dir="$1"
    local label="$2"

    grep "${label}" "${case_dir}/${SOLVER_LOGFILE}" 2>/dev/null \
        | tail -n 1 \
        | awk '{print $NF}' \
        || true
}

check_less_than() {
    local case_name="$1"
    local label="$2"
    local value="$3"
    local tolerance="$4"

    if [[ -z "${value}" ]]; then
        echo "FAIL: ${case_name}: could not extract ${label}"
        failures=$((failures + 1))
    elif awk "BEGIN {exit !(${value} < ${tolerance})}"; then
        printf "PASS: %s: %s = %.6g\n" "${case_name}" "${label}" "${value}"
    else
        printf "FAIL: %s: %s = %.6g\n" "${case_name}" "${label}" "${value}"
        failures=$((failures + 1))
    fi
}

compare_written_fields() {
    python3 - "$1" "$2" << 'PYEOF'
import re
import sys

number = re.compile(
    r"(?<![A-Za-z_])[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
)

try:
    texts = [open(path).read() for path in sys.argv[1:]]
    structure = [number.sub("<number>", text) for text in texts]
    if structure[0] != structure[1]:
        raise ValueError("field structures differ")

    values = [
        [float(match.group()) for match in number.finditer(text)]
        for text in texts
    ]
    if not values[0] or len(values[0]) != len(values[1]):
        raise ValueError("field numeric data are missing or differ in size")

    print(max(abs(a - b) for a, b in zip(*values)))
except (OSError, ValueError) as error:
    print(error, file=sys.stderr)
    sys.exit(1)
PYEOF
}

# The arms whose solidProperties carries useMechanicalConstitutiveLawManager.
# Every other arm must run the legacy path
FRAMEWORK_APPROACHES=(
    segregatedManager
    petscSnesPressureManager
    highOrderFourthOrder
    pressureDisplacementCompressibleManager
    pressureDisplacementIncompressibleManager
)

is_framework_approach() {
    local a="$1"
    local f
    for f in "${FRAMEWORK_APPROACHES[@]}"; do
        [[ "${a}" == "${f}" ]] && return 0
    done
    return 1
}

# An arm that silently ran the other path still satisfies every tolerance
# below, because both paths solve the same problem correctly. Without this the
# framework arms prove nothing: losing the switch from the dictionary, or
# linking the wrong file, would read as a pass
check_took_its_path() {
    local approach="$1"
    local log="${2}/${SOLVER_LOGFILE}"

    if [[ ! -f "${log}" ]]; then
        echo "FAIL: ${approach}: no solver log to check the path taken"
        failures=$((failures + 1))
        return
    fi

    if is_framework_approach "${approach}"; then
        if ! grep -q "Selecting mechanical constitutive law" "${log}"; then
            echo "FAIL: ${approach}: framework arm did not use the framework"
            failures=$((failures + 1))
        fi
    else
        if grep -q "Selecting mechanical constitutive law" "${log}"; then
            echo "FAIL: ${approach}: legacy arm used the framework"
            failures=$((failures + 1))
        fi
    fi
}

failures=0

for approach in "${APPROACHES[@]}"; do
    case_dir="${REGRESSION_ROOT}/${approach}"
    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: ${approach}"
        continue
    fi
    check_took_its_path "${approach}" "${case_dir}"
    check_less_than \
        "${approach}" "DDifference LInf" \
        "$(extract_disp_linf "${case_dir}" "DDifference")" \
        "${DISP_TOL}"
    check_less_than \
        "${approach}" "pointDDifference LInf" \
        "$(extract_disp_linf "${case_dir}" "pointDDifference")" \
        "${POINT_DISP_TOL}"
    check_less_than \
        "${approach}" "stress component-0 LInf" \
        "$(extract_stress_linf_comp0 "${case_dir}")" \
        "${STRESS_TOL}"
done

# ------------------------------------------------------------
# Each framework arm against the legacy arm it mirrors
# ------------------------------------------------------------
# Checking both arms against the analytical tolerances separately does not
# compare them with each other: both paths solve this problem well within
# tolerance, so both would pass even if they disagreed. For isotropic linear
# elasticity a deviatoric projection and the declared volumetric split are the
# same operation, so the two arms must agree to the precision written in the
# fields, not merely satisfy the analytical tolerances independently.
#
# highOrderFourthOrder has no legacy twin here - the other high-order arms are
# different discretisations rather than the same one on the legacy path - so it
# is covered by its tolerances and its path assertion only
FRAMEWORK_PAIRS=(
    "segregated segregatedManager"
    "petscSnesPressure petscSnesPressureManager"
)

for pair in "${FRAMEWORK_PAIRS[@]}"; do
    IFS=' ' read -r legacy_arm framework_arm <<< "${pair}"
    legacy_case="${REGRESSION_ROOT}/${legacy_arm}"
    framework_case="${REGRESSION_ROOT}/${framework_arm}"

    if solids4Foam::regressionCaseSkipped "${legacy_case}/${ALLRUN_LOGFILE}" \
        || solids4Foam::regressionCaseSkipped \
            "${framework_case}/${ALLRUN_LOGFILE}"
    then
        echo "SKIP: ${framework_arm} against ${legacy_arm}"
        continue
    fi

    t=$(solids4Foam::latestTime "${framework_case}")

    if [[ -z "${t}" || ! -f "${framework_case}/${t}/D" \
        || ! -f "${legacy_case}/${t}/D" ]]
    then
        echo "FAIL: ${framework_arm}: no D field to compare with ${legacy_arm}"
        failures=$((failures + 1))
        continue
    fi

    field_diff=""
    if field_diff=$(compare_written_fields \
        "${legacy_case}/${t}/D" "${framework_case}/${t}/D") \
      && awk "BEGIN {exit !(${field_diff} <= ${FRAMEWORK_FIELD_ABS_TOL})}"
    then
        printf "PASS: %s and %s agree to write precision (max |delta| = %.3g)\n" \
            "${framework_arm}" "${legacy_arm}" "${field_diff}"
    else
        printf "FAIL: %s and %s fields differ (max |delta| = %s)\n" \
            "${framework_arm}" "${legacy_arm}" "${field_diff:-unavailable}"
        failures=$((failures + 1))
    fi
done

for case_args in "${PRESSURE_DISPLACEMENT_CASES[@]}"; do
    IFS=' ' read -r approach mesh <<< "${case_args}"
    case_name="${approach}-${mesh}"
    case_dir="${REGRESSION_ROOT}/${case_name}"

    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: ${case_name}"
        continue
    fi

    check_took_its_path "${approach}" "${case_dir}"

    if [[ "${approach}" == *Manager ]]; then
        twin_dir="${REGRESSION_ROOT}/${approach%Manager}-${mesh}"
        for label in "DError, max" "pErr, max"; do
            twin_value="$(extract_log_value "${twin_dir}" "${label}")"
            value="$(extract_log_value "${case_dir}" "${label}")"
            if [[ -z "${twin_value}" || -z "${value}" ]]; then
                echo "FAIL: ${case_name}: could not compare ${label} with its twin"
                failures=$((failures + 1))
            elif awk "BEGIN {d = ${value} - ${twin_value}; \
                exit !(${twin_value} > 0 \
                    && d*d <= (${PD_FRAMEWORK_REL_TOL}*${twin_value})^2)}"
            then
                echo "PASS: ${case_name}: ${label} matches the legacy twin"
            else
                echo "FAIL: ${case_name}: ${label} = ${value}, legacy twin ${twin_value}"
                failures=$((failures + 1))
            fi
        done
    fi

    check_less_than \
        "${case_name}" "DError, max" \
        "$(extract_log_value "${case_dir}" "DError, max")" "${PD_DISP_TOL}"
    check_less_than \
        "${case_name}" "pointDError, max" \
        "$(extract_log_value "${case_dir}" "pointDError, max")" \
        "${PD_POINT_DISP_TOL}"
    check_less_than \
        "${case_name}" "sigmaXXErr, max" \
        "$(extract_log_value "${case_dir}" "sigmaXXErr, max")" \
        "${PD_SIGMA_TOL}"
    check_less_than \
        "${case_name}" "sigmaXYErr, max" \
        "$(extract_log_value "${case_dir}" "sigmaXYErr, max")" \
        "${PD_SIGMA_TOL}"
    check_less_than \
        "${case_name}" "sigmaYYErr, max" \
        "$(extract_log_value "${case_dir}" "sigmaYYErr, max")" \
        "${PD_SIGMA_TOL}"
    check_less_than \
        "${case_name}" "pErr, max" \
        "$(extract_log_value "${case_dir}" "pErr, max")" "${PD_P_TOL}"
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
