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
# suctionCaission regression test
#
# Pullout of a suction caisson: poroMechanicalLaw over
# linearElasticMohrCoulombPlastic, on poroLinearGeometry.
#
# The bands are wide, and deliberately so. This case does not
# converge tightly - it reaches the corrector limit on every
# time step, as it did with the removed legacy mechanicalModel
# - so the bands say the caisson yielded and suction developed,
# not that a particular number came back.
# ============================================================

EPS_MIN=0.85
EPS_MAX=1.25
SIGMA_MIN=1.5e6
SIGMA_MAX=2.2e6

# The comparison with the removed legacy mechanicalModel runs to t = 1 rather
# than ten. The full case takes about nine minutes, and the plastic return
# mapping and the effective stress the composite carries are both exercised
# within the first step
COMPARISON_END_TIME=1

# The legacy model's D and porePressure at t = 1, as the max and mean component
# magnitude of each field written to fourteen figures, from the last commit
# that had it (mcl-stage8-coverage, c3a92b3d), per fork. The framework
# reproduced both fields exactly, in every one of those figures. These are
# recorded numbers, though, and another compiler, CPU or MPI build moves an
# iterative solution by round-off at the solver tolerance: CI measures up to
# 3e-8 relative against values recorded on macOS. The tolerance, 1e-6 of the
# largest value, allows for that, and is ten times below the 1e-5 that a
# 0.001% change in a material constant makes
case "$(solids4Foam::foamFlavour)" in
    com)
        LEGACY_D_MAX=0.039970270753559
        LEGACY_D_MEAN=0.00324034872160463
        LEGACY_P_MAX=94941.676957032
        LEGACY_P_MEAN=19655.7130753187
        ;;
    org)
        LEGACY_D_MAX=0.039970258883115
        LEGACY_D_MEAN=0.00324070141058448
        LEGACY_P_MAX=94957.392747739
        LEGACY_P_MEAN=19661.8147230215
        ;;
    foamextend)
        LEGACY_D_MAX=0.039984257331651
        LEGACY_D_MEAN=0.00342262334429839
        LEGACY_P_MAX=66272.071942348
        LEGACY_P_MEAN=17304.5184317787
        ;;
esac
LEGACY_REL_TOL=1e-6

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "suctionCaission regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "Plus the comparison with the legacy model, to t=${COMPARISON_END_TIME}"
echo "============================================================"
echo

prepare_case() {
    local dir="$1"

    rm -rf "${dir}"
    mkdir -p "${dir}"

    local item base_item
    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${dir}/"
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
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped here"
    exit 0
fi

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

# A field against the removed legacy model's, through the norms above. Both
# differences are bounded by the largest pointwise difference, so a field that
# agrees with the legacy one to tol times its largest value passes, and one
# that does not is caught by at least one of the two in all but contrived cases
check_field_against_legacy() {
    local label="$1"
    local file="$2"
    local legacy_max="$3"
    local legacy_mean="$4"
    local tol="$5"
    local norms field_max field_mean

    if [[ ! -f "${file}" ]] || ! norms=$(internal_field_norms "${file}"); then
        echo "FAIL: ${label}: no field to compare with the legacy model"
        return 1
    fi

    read -r field_max field_mean <<< "${norms}"

    if awk "BEGIN {
            a = ${field_max} - ${legacy_max}; if (a < 0) a = -a
            b = ${field_mean} - ${legacy_mean}; if (b < 0) b = -b
            exit !(${field_max} > 0 && a <= ${tol}*${legacy_max} \
                && b <= ${tol}*${legacy_max})
        }"
    then
        printf "PASS: %s matches the legacy model: max %.15g (%.15g), mean %.15g (%.15g)\n" \
            "${label}" "${field_max}" "${legacy_max}" "${field_mean}" "${legacy_mean}"
        return 0
    fi

    printf "FAIL: %s differs from the legacy model: max %.15g (%.15g), mean %.15g (%.15g), tolerance %s\n" \
        "${label}" "${field_max}" "${legacy_max}" "${field_mean}" "${legacy_mean}" "${tol}"
    return 1
}

# Run the case to the comparison time and hold it to the legacy answer there
run_legacy_comparison() {
    local dir="${REGRESSION_ROOT}/comparison"

    prepare_case "${dir}"

    sed -i "s|^endTime .*|endTime         ${COMPARISON_END_TIME};|" \
        "${dir}/system/controlDict"

    # Enough digits that the comparison is about the solution and not about
    # the last figure written
    if grep -q "^writePrecision" "${dir}/system/controlDict"; then
        sed -i 's|^writePrecision.*|writePrecision  14;|' \
            "${dir}/system/controlDict"
    else
        echo "writePrecision  14;" >> "${dir}/system/controlDict"
    fi

    ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
        echo "FAIL: the comparison could not run ${dir}"
        return 1
    }

    if ! grep -q "Selecting mechanical constitutive law" \
        "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the comparison constructed no mechanical constitutive law"
        return 1
    fi

    local t
    t=$(solids4Foam::latestTime "${dir}")

    if [[ -z "${t}" ]] \
        || ! awk "BEGIN {exit !((${t} - ${COMPARISON_END_TIME})^2 <= 1e-20)}"
    then
        echo "FAIL: the comparison stopped at '${t}', not at ${COMPARISON_END_TIME}"
        return 1
    fi

    local failed=0

    check_field_against_legacy "D at t = ${t}" "${dir}/${t}/D" \
        "${LEGACY_D_MAX}" "${LEGACY_D_MEAN}" "${LEGACY_REL_TOL}" \
        || failed=1

    check_field_against_legacy "porePressure at t = ${t}" \
        "${dir}/${t}/porePressure" \
        "${LEGACY_P_MAX}" "${LEGACY_P_MEAN}" "${LEGACY_REL_TOL}" \
        || failed=1

    return "${failed}"
}

epsilon=$(grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null \
    | tail -n 1 | awk '{print $NF}' || true)
sigma=$(grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
    2>/dev/null | tail -n 1 | awk '{print $NF}' || true)

if [[ -z "${epsilon}" || -z "${sigma}" ]]; then
    echo "FAIL: could not extract the regression quantities"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${epsilon} >= ${EPS_MIN} && ${epsilon} <= ${EPS_MAX})}"
then
    printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
else
    printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${sigma} <= ${SIGMA_MAX})}"
then
    printf "PASS: Max sigmaEq = %.6g\n" "${sigma}"
else
    printf "FAIL: Max sigmaEq = %.6g\n" "${sigma}"
    failures=$((failures + 1))
fi

if [ "$CHECK_ONLY" = false ] && ! run_legacy_comparison; then
    failures=$((failures + 1))
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
