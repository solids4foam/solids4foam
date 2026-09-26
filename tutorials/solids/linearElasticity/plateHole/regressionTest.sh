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
# (segregated/petscSnes/petscSnesPressure/high-order variants) and
# pressure-displacement
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

PD_DISP_TOL=3.0e-4
PD_POINT_DISP_TOL=3.0e-4
PD_SIGMA_TOL=2.0e5
PD_P_TOL=9.0e4

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
PARALLEL_N_PROCS=2

APPROACHES=(
    segregated
    petscSnes
    petscSnesPressure
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
)

# ------------------------------------------------------------
# Against the removed legacy mechanicalModel
# ------------------------------------------------------------
# From the last commit that had it (mcl-stage8-coverage, c3a92b3d).
#
# segregated and petscSnesPressure: the final D, as the max and mean component
# magnitude of the field as written, per fork. For isotropic linear elasticity
# a deviatoric projection and the declared volumetric split are the same
# operation, so the framework reproduced the legacy D fields to the precision
# written, 2e-12 at most. D is written to six figures, though, so a round-off
# difference on another compiler, CPU or MPI build can move the last one; the
# tolerance, 2e-5 of the largest value, is two units in that figure. Both norms
# are within the largest pointwise difference, so the bound carries over
declare -A LEGACY_D_MAX=()
declare -A LEGACY_D_MEAN=()
case "$(solids4Foam::foamFlavour)" in
    com|org)
        LEGACY_D_MAX[segregated]=1.05063e-05
        LEGACY_D_MEAN[segregated]=2.56367713626666e-06
        LEGACY_D_MAX[petscSnesPressure]=1.05169e-05
        LEGACY_D_MEAN[petscSnesPressure]=2.56525614836667e-06
        ;;
    foamextend)
        LEGACY_D_MAX[segregated]=1.05101e-05
        LEGACY_D_MEAN[segregated]=2.56399595763334e-06
        ;;
esac
LEGACY_D_REL_TOL=2e-5

# The pressure-displacement cases, on foam-extend 4.1, where
# coupledPressureDisplacementSolid runs: DError and pErr maxima, as logged. The
# model runs in linear mode here, where the law is not evaluated and the
# stiffness is the same shear modulus on both paths, so the framework matched
# them to the 1e-6 relative it was held to. These are the legacy values CI
# logged (foam-extend-4.1-PETSc image, mcl-stage8-coverage 652cdb62). A macOS
# foam-extend 4.1 build gave 0.0001275 and 38740.5 for the compressible case,
# up to 40% apart, so this solver's answer here depends on the build, not only
# on round-off, and the values hold for the CI build. The log gives six
# figures; the tolerance, 1e-5 relative, is a few units in the last one
declare -A LEGACY_PD=(
    # case : DError max, pErr max
    [pressureDisplacementCompressible-coarse]="7.748e-05 22215.6"
    [pressureDisplacementIncompressible-coarse]="8.66026e-05 31612.4"
)
LEGACY_PD_REL_TOL=1e-5

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

# The largest magnitude of any component of a field's internal values, and the
# mean magnitude, as "max<TAB>mean"
internal_field_norms() {
    python3 - "$1" << 'PYEOF'
import re
import sys

number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
text = open(sys.argv[1]).read()

uniform = re.search(
    r"\binternalField\s+uniform\s+(\([^)]*\)|" + number + r")\s*;", text
)
if uniform:
    body = uniform.group(1)
else:
    field = re.search(
        r"\binternalField\s+nonuniform\s+List<\w+>\s+\d+\s*\((.*?)\n\)\s*;",
        text,
        re.DOTALL,
    )
    if not field:
        sys.exit(f"cannot parse internalField in {sys.argv[1]}")
    body = field.group(1)

values = [abs(float(x)) for x in re.findall(number, body)]
if not values:
    sys.exit(f"empty internalField in {sys.argv[1]}")
print(f"{max(values):.15g}\t{sum(values)/len(values):.15g}")
PYEOF
}


# Every arm takes its material from the mechanicalConstitutiveLaw framework
check_used_framework() {
    local approach="$1"
    local log="${2}/${SOLVER_LOGFILE}"

    if [[ ! -f "${log}" ]]; then
        echo "FAIL: ${approach}: no solver log to check"
        failures=$((failures + 1))
    elif ! grep -q "Selecting mechanical constitutive law" "${log}"; then
        echo "FAIL: ${approach}: constructed no mechanical constitutive law"
        failures=$((failures + 1))
    fi
}

failures=0

for approach in "${APPROACHES[@]}"; do
    case_dir="${REGRESSION_ROOT}/${approach}"
    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: ${approach}"
        continue
    fi
    check_used_framework "${approach}" "${case_dir}"
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
# segregated and petscSnesPressure against the legacy D fields
# ------------------------------------------------------------
# Checking against the analytical tolerances does not do this: both models
# solve this problem well within them, so a framework that disagreed with the
# legacy model would still pass them.
#
# highOrderFourthOrder has no legacy counterpart - it asks for the full
# material tangent, which the legacy model did not have - so it is covered by
# its tolerances only
for approach in segregated petscSnesPressure; do
    case_dir="${REGRESSION_ROOT}/${approach}"

    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: ${approach} against the legacy model"
        continue
    fi

    t=$(solids4Foam::latestTime "${case_dir}")
    legacy_max="${LEGACY_D_MAX[${approach}]:-}"
    legacy_mean="${LEGACY_D_MEAN[${approach}]:-}"

    if [[ -z "${legacy_max}" ]]; then
        echo "FAIL: ${approach}: no legacy answer recorded for this fork"
        failures=$((failures + 1))
        continue
    fi

    if [[ -z "${t}" || ! -f "${case_dir}/${t}/D" ]] \
        || ! norms=$(internal_field_norms "${case_dir}/${t}/D")
    then
        echo "FAIL: ${approach}: no D field to compare with the legacy model"
        failures=$((failures + 1))
        continue
    fi

    read -r field_max field_mean <<< "${norms}"

    if awk "BEGIN {
            a = ${field_max} - ${legacy_max}; if (a < 0) a = -a
            b = ${field_mean} - ${legacy_mean}; if (b < 0) b = -b
            tol = ${LEGACY_D_REL_TOL}*${legacy_max}
            exit !(${field_max} > 0 && a <= tol && b <= tol)
        }"
    then
        printf "PASS: %s D matches the legacy model (max %.10g, mean %.10g)\n" \
            "${approach}" "${field_max}" "${field_mean}"
    else
        printf "FAIL: %s D differs from the legacy model: max %.10g (%.10g), mean %.10g (%.10g)\n" \
            "${approach}" "${field_max}" "${legacy_max}" "${field_mean}" \
            "${legacy_mean}"
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

    if ! grep -q "taking the stiffness from the mechanicalConstitutiveLaw" \
        "${case_dir}/${SOLVER_LOGFILE}" 2>/dev/null
    then
        echo "FAIL: ${case_name}: did not take its stiffness from the framework"
        failures=$((failures + 1))
    fi

    if [[ -n "${LEGACY_PD[${case_name}]:-}" ]]; then
        IFS=' ' read -r legacy_derror legacy_perr \
            <<< "${LEGACY_PD[${case_name}]}"
        for pair in "DError, max:${legacy_derror}" "pErr, max:${legacy_perr}"
        do
            label="${pair%%:*}"
            legacy_value="${pair##*:}"
            value="$(extract_log_value "${case_dir}" "${label}")"
            if [[ -z "${value}" ]]; then
                echo "FAIL: ${case_name}: could not compare ${label} with the legacy model"
                failures=$((failures + 1))
            elif awk "BEGIN {d = ${value} - ${legacy_value}; \
                exit !(${legacy_value} > 0 \
                    && d*d <= (${LEGACY_PD_REL_TOL}*${legacy_value})^2)}"
            then
                echo "PASS: ${case_name}: ${label} matches the legacy model"
            else
                echo "FAIL: ${case_name}: ${label} = ${value}, legacy model ${legacy_value}"
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
