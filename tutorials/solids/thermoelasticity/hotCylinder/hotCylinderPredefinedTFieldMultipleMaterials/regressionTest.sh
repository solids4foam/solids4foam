#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# hotCylinderPredefinedTFieldMultipleMaterials regression test
#
# Two materials, steel and aluminium, whose temperature is not solved for but
# read at each time step from the separate case hotCylinderTemperatureField
# (TcaseDirectory). Two arms, each held to the answer of the removed legacy
# mechanicalModel:
#
#   main              the tutorial as it is
#   single            steel alone
#
# The single-material arm isolates the reading of T. With one material the
# framework discretised the problem as the legacy model did, so any difference
# beyond the solution tolerance would be the temperature it was given. The
# two-material arm is then the multi-material path with that input.
# ============================================================

# Bounds on the main arm's final extrema. The temperature is uniform in
# each time directory, so the stress is the mismatch of the two materials'
# expansion. Measured on the legacy model: max epsilonEq 1.47489e-3 and max
# sigmaEq 1.8164e8 on OpenFOAM.com v2512 and OpenFOAM.org 9, 1.47582e-3 and
# 1.81641e8 on foam-extend 4.1
EPS_MIN=1.4e-3
EPS_MAX=1.55e-3
SIGMA_MIN=1.75e8
SIGMA_MAX=1.9e8

# The final D of the removed legacy mechanicalModel, as the max and mean
# component magnitude of the field written to fourteen figures, for each arm,
# from the last commit that had it (mcl-stage8-coverage, c3a92b3d), per fork.
case "$(solids4Foam::foamFlavour)" in
    com)
        LEGACY_D_MAX=0.00098799181372038
        LEGACY_D_MEAN=0.000340730093870825
        LEGACY_SINGLE_D_MAX=0.0007175233920732
        LEGACY_SINGLE_D_MEAN=0.000264826849529889
        ;;
    org)
        LEGACY_D_MAX=0.00098799184046576
        LEGACY_D_MEAN=0.000340730094882777
        LEGACY_SINGLE_D_MAX=0.00071752343529793
        LEGACY_SINGLE_D_MEAN=0.000264826849186313
        ;;
    foamextend)
        LEGACY_D_MAX=0.00098782516201107
        LEGACY_D_MEAN=0.000340707432130164
        LEGACY_SINGLE_D_MAX=0.00071750897253875
        LEGACY_SINGLE_D_MEAN=0.000264827217209345
        ;;
esac

# Two materials. These are different discretisations of the interface - per
# material sub-meshes on the legacy path, the material-aware leastSquaresS4f
# gradient on one mesh on the framework - as on layeredPipe and punch, so they
# agree only to the discretisation. Measured: 1.23e-3 of the largest
# displacement on OpenFOAM.com v2512 and OpenFOAM.org 9, 1.35e-3 on
# foam-extend 4.1. This is the coarse check: a 1% change in steel's alpha
# took the framework to 6.2e-3 and failed, a 0.1% change (1.7e-3) did not, and
# holding T at its second-step value took it to 0.12. The single-material arm
# below is the fine one
FRAMEWORK_D_REL_TOL=3e-3

# One material. The same discretisation, so the two agreed to the solution
# tolerance (1e-6): measured 2.7e-8 on OpenFOAM.com v2512, 3.1e-8 on
# foam-extend 4.1 and 8.2e-8 on OpenFOAM.org 9. The threshold is well above
# that and ten times below the 1.0e-5 that a 0.001% change in alpha made
SINGLE_D_REL_TOL=1e-6

# The last time step, and the last directory of the temperature case
COMPARISON_END_TIME=4

# A thermal stress of at least this much must be present at the end in every
# arm, so that the comparisons are not of stress-free solutions
SIGMA_NONZERO_MIN=1e7

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
FRAMEWORK_MARK="Selecting mechanical constitutive law"
FRAMEWORK_T_READ="Reading T for the mechanical constitutive laws from"

echo "============================================================"
echo "hotCylinderPredefinedTFieldMultipleMaterials regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "D against the legacy model, relative to its largest value: two"\
" materials < ${FRAMEWORK_D_REL_TOL}, one material < ${SINGLE_D_REL_TOL}"
echo "============================================================"
echo

copy_case() {
    local dest="$1"

    rm -rf "${dest}"
    mkdir -p "${dest}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${dest}/"
    done

    # Enough digits that the comparisons measure the solutions rather than
    # the last digit written
    sed -i 's|^writePrecision.*|writePrecision  14;|' \
        "${dest}/system/controlDict"
}

# Keep steel only. One law covers the whole mesh, so no cellZone is used
use_single_material() {
    local dir="$1"

    python3 - "${dir}/constant/mechanicalProperties" << 'PYEOF'
import sys
path = sys.argv[1]
text = open(path).read()
start = text.index("    aluminium")
end = text.index(");", start)
open(path, "w").write(text[:start] + text[end:])
PYEOF

    if grep -q aluminium "${dir}/constant/mechanicalProperties"; then
        echo "FAIL: could not remove the second material from ${dir}"
        return 1
    fi

    # One material needs no material-aware gradient, and the legacy single
    # material answer this arm is held to was computed with the gradient the
    # tutorial used before it needed one. Using it here keeps the comparison
    # about the temperature the law was given, not about the gradient
    sed -i \
        's|^\( *default *\)leastSquaresS4f;|\1pointCellsLeastSquares;|' \
        "${dir}/system/fvSchemes"

    if ! grep -q "pointCellsLeastSquares" "${dir}/system/fvSchemes"; then
        echo "FAIL: could not set the single-material gradient in ${dir}"
        return 1
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

extract_last() {
    grep "$1" "$2/${SOLVER_LOGFILE}" | tail -n 1 | awk '{print $NF}'
}

# Checks common to every arm: it completed, took its material and its
# temperature from the framework, and ended at the last time with a real
# thermal stress
check_arm() {
    local dir="$1"

    if ! grep -q "^End" "${dir}/${SOLVER_LOGFILE}" 2>/dev/null \
      || grep -qE "FOAM FATAL|did not converge" "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${dir##*/} did not complete and converge"
        return 1
    fi

    if ! grep -q "${FRAMEWORK_MARK}" "${dir}/${SOLVER_LOGFILE}"; then
        echo "FAIL: ${dir##*/} constructed no mechanical constitutive law"
        return 1
    fi

    # The input must have come from the temperature case, and at the last time
    # too, not merely once at the start
    if ! grep "${FRAMEWORK_T_READ}" "${dir}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | grep -q "hotCylinderTemperatureField/${COMPARISON_END_TIME}\""
    then
        echo "FAIL: ${dir##*/} did not read T from" \
            "hotCylinderTemperatureField/${COMPARISON_END_TIME}"
        return 1
    fi

    local t
    t=$(solids4Foam::latestTime "${dir}")

    if ! awk "BEGIN {exit !((${t:-0} - ${COMPARISON_END_TIME})^2 <= 1e-20)}"
    then
        echo "FAIL: ${dir##*/} stopped at '${t}';" \
            "expected ${COMPARISON_END_TIME}"
        return 1
    fi

    local sigma
    sigma=$(extract_last "Max sigmaEq (von Mises stress)" "${dir}")

    if [[ -z "${sigma}" ]] \
      || ! awk "BEGIN {exit !(${sigma} >= ${SIGMA_NONZERO_MIN})}"
    then
        echo "FAIL: ${dir##*/} has no thermal stress (max sigmaEq '${sigma}')"
        return 1
    fi
}

# ------------------------------------------------------------
# Clean & run the arms
# ------------------------------------------------------------

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

SINGLE_DIR="${REGRESSION_ROOT}/single"

if [ "$CHECK_ONLY" = false ]; then
    for dir in "${CASE_DIR}" "${SINGLE_DIR}"; do
        copy_case "${dir}"
    done

    use_single_material "${SINGLE_DIR}"

    for dir in "${CASE_DIR}" "${SINGLE_DIR}"; do
        ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check_arm "${CASE_DIR}" || failures=$((failures + 1))
check_arm "${SINGLE_DIR}" || failures=$((failures + 1))

epsilon=$(extract_last "Max epsilonEq" "${CASE_DIR}" || true)
sigma=$(extract_last "Max sigmaEq (von Mises stress)" "${CASE_DIR}" || true)

if [[ -z "${epsilon}" || -z "${sigma}" ]]; then
    echo "FAIL: Could not extract one or more regression quantities"
    failures=$((failures + 1))
else
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
fi

check_field_against_legacy "two materials, D" \
    "${CASE_DIR}/${COMPARISON_END_TIME}/D" \
    "${LEGACY_D_MAX}" "${LEGACY_D_MEAN}" "${FRAMEWORK_D_REL_TOL}" \
    || failures=$((failures + 1))

check_field_against_legacy "one material, D" \
    "${SINGLE_DIR}/${COMPARISON_END_TIME}/D" \
    "${LEGACY_SINGLE_D_MAX}" "${LEGACY_SINGLE_D_MEAN}" "${SINGLE_D_REL_TOL}" \
    || failures=$((failures + 1))

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
    for dir in "${CASE_DIR}" "${SINGLE_DIR}"; do
        ( cd "${dir}" && ./Allclean > /dev/null 2>&1 ) || true
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
