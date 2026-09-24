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
# (TcaseDirectory). Four arms:
#
#   main              the tutorial as it is, on the legacy mechanicalModel
#   framework         the same on the mechanicalConstitutiveLaw framework
#   single            steel alone, on the legacy mechanicalModel
#   singleFramework   steel alone, on the framework
#
# The single-material pair isolates the reading of T. With one material both
# paths discretise the problem identically, so any difference beyond the
# solution tolerance would be the temperature they were given. The
# two-material pair is then the multi-material path with that input.
# ============================================================

# Bounds on the legacy arm's final extrema. The temperature is uniform in
# each time directory, so the stress is the mismatch of the two materials'
# expansion. Measured: max epsilonEq 1.47489e-3 and max sigmaEq 1.8164e8 on
# OpenFOAM.com v2512 and OpenFOAM.org 9, 1.47582e-3 and 1.81641e8 on
# foam-extend 4.1
EPS_MIN=1.4e-3
EPS_MAX=1.55e-3
SIGMA_MIN=1.75e8
SIGMA_MAX=1.9e8

# The framework arm against the legacy arm, two materials. These are
# different discretisations of the interface - per-material sub-meshes on the
# legacy path, the material-aware leastSquaresS4f gradient on one mesh on the
# framework path - as on layeredPipe and punch, so they agree only to the
# discretisation. Measured: 1.23e-3 of the largest displacement on
# OpenFOAM.com v2512 and OpenFOAM.org 9, 1.35e-3 on foam-extend 4.1. This is
# the coarse check: a 1% change in steel's alpha on the framework arm takes it
# to 6.2e-3 and fails, a 0.1% change (1.7e-3) does not, and holding T at its
# second-step value takes it to 0.12. The single-material pair below is the
# fine one
FRAMEWORK_D_REL_TOL=3e-3

# The framework arm against the legacy arm, one material. The same
# discretisation, so they agree to the solution tolerance (1e-6): measured
# 2.7e-8 on OpenFOAM.com v2512, 3.1e-8 on foam-extend 4.1 and 8.2e-8 on
# OpenFOAM.org 9. The threshold is well above that and ten times below the
# 1.0e-5 that a 0.001% change in alpha on the framework arm makes
SINGLE_D_REL_TOL=1e-6

# The last time step, and the last directory of the temperature case
COMPARISON_END_TIME=4

# A thermal stress of at least this much must be present at the end in every
# arm, so that the comparisons are not of four stress-free solutions
SIGMA_NONZERO_MIN=1e7

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
FRAMEWORK_MARK="Selecting mechanical constitutive law"
FRAMEWORK_T_READ="Reading T for the mechanical constitutive laws from"

echo "============================================================"
echo "hotCylinderPredefinedTFieldMultipleMaterials regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "Framework vs legacy relative D difference: two materials < "\
"${FRAMEWORK_D_REL_TOL}, one material < ${SINGLE_D_REL_TOL}"
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

# The switch goes in the solid model's coeffs sub-dictionary
use_framework() {
    local dir="$1"

    sed -i \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${dir}/constant/solidProperties"
}

# More than one material on the framework needs the material-aware gradient.
# Only then: the single-material framework arm keeps the legacy arm's
# gradient, so that the pair differ in nothing but the path
use_material_aware_gradient() {
    local dir="$1"

    sed -i \
        's|^\( *default *\)pointCellsLeastSquares;|\1leastSquaresS4f;|' \
        "${dir}/system/fvSchemes"
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
}

compare_internal_vector_fields() {
    python3 - "$1" "$2" << 'PYEOF'
import re
import sys

number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

def read_internal(path):
    text = open(path).read()
    field = re.search(
        r"\binternalField\s+nonuniform\s+List<vector>\s+\d+\s*\((.*?)\)\s*;",
        text,
        re.DOTALL,
    )
    if not field:
        raise ValueError(f"cannot parse a nonuniform internalField in {path}")
    values = re.findall(
        rf"\(({number})\s+({number})\s+({number})\)", field.group(1)
    )
    if not values:
        raise ValueError(f"empty internalField in {path}")
    return [tuple(map(float, value)) for value in values]

try:
    a = read_internal(sys.argv[1])
    b = read_internal(sys.argv[2])
    if len(a) != len(b):
        raise ValueError("different internalField sizes")
    max_diff = max(abs(x - y) for av, bv in zip(a, b) for x, y in zip(av, bv))
    max_value = max(abs(x) for av in a for x in av)
    if max_value == 0.0:
        raise ValueError("the reference displacement is zero")
    print(f"{max_diff / max_value:.12g}")
except Exception as exc:
    print(exc, file=sys.stderr)
    sys.exit(1)
PYEOF
}

extract_last() {
    grep "$1" "$2/${SOLVER_LOGFILE}" | tail -n 1 | awk '{print $NF}'
}

# Checks common to every arm: it completed, took the path it was set up for,
# and ended at the last time with a real thermal stress
check_arm() {
    local dir="$1"
    local expect_framework="$2"

    if ! grep -q "^End" "${dir}/${SOLVER_LOGFILE}" 2>/dev/null \
      || grep -qE "FOAM FATAL|did not converge" "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${dir##*/} did not complete and converge"
        return 1
    fi

    if [[ "${expect_framework}" == yes ]]; then
        if ! grep -q "${FRAMEWORK_MARK}" "${dir}/${SOLVER_LOGFILE}"; then
            echo "FAIL: ${dir##*/} did not use the framework"
            return 1
        fi

        # The input must have come from the temperature case, and at the
        # last time too, not merely once at the start
        if ! grep "${FRAMEWORK_T_READ}" "${dir}/${SOLVER_LOGFILE}" \
            | tail -n 1 \
            | grep -q "hotCylinderTemperatureField/${COMPARISON_END_TIME}"
        then
            echo "FAIL: ${dir##*/} did not read T from" \
                "hotCylinderTemperatureField/${COMPARISON_END_TIME}"
            return 1
        fi
    elif grep -q "${FRAMEWORK_MARK}" "${dir}/${SOLVER_LOGFILE}"; then
        echo "FAIL: ${dir##*/} used the framework"
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

compare_arms() {
    local reference="$1"
    local candidate="$2"
    local tol="$3"
    local label="$4"
    local f="${COMPARISON_END_TIME}/D"

    if [[ ! -f "${reference}/${f}" || ! -f "${candidate}/${f}" ]]; then
        echo "FAIL: ${label}: no ${f} to compare"
        return 1
    fi

    local rel
    if ! rel=$(compare_internal_vector_fields \
        "${reference}/${f}" "${candidate}/${f}")
    then
        echo "FAIL: ${label}: could not compare the D internal fields"
        return 1
    fi

    if awk "BEGIN {exit !(${rel} < ${tol})}"; then
        printf "PASS: %s: relative D diff = %.4g (< %s)\n" \
            "${label}" "${rel}" "${tol}"
        return 0
    fi

    printf "FAIL: %s: relative D diff = %.4g (>= %s)\n" \
        "${label}" "${rel}" "${tol}"
    return 1
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

FRAMEWORK_DIR="${REGRESSION_ROOT}/framework"
SINGLE_DIR="${REGRESSION_ROOT}/single"
SINGLE_FRAMEWORK_DIR="${REGRESSION_ROOT}/singleFramework"

if [ "$CHECK_ONLY" = false ]; then
    for dir in "${CASE_DIR}" "${FRAMEWORK_DIR}" "${SINGLE_DIR}" \
        "${SINGLE_FRAMEWORK_DIR}"
    do
        copy_case "${dir}"
    done

    use_framework "${FRAMEWORK_DIR}"
    use_material_aware_gradient "${FRAMEWORK_DIR}"
    use_single_material "${SINGLE_DIR}"
    use_single_material "${SINGLE_FRAMEWORK_DIR}"
    use_framework "${SINGLE_FRAMEWORK_DIR}"

    for dir in "${CASE_DIR}" "${FRAMEWORK_DIR}" "${SINGLE_DIR}" \
        "${SINGLE_FRAMEWORK_DIR}"
    do
        ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check_arm "${CASE_DIR}" no || failures=$((failures + 1))
check_arm "${FRAMEWORK_DIR}" yes || failures=$((failures + 1))
check_arm "${SINGLE_DIR}" no || failures=$((failures + 1))
check_arm "${SINGLE_FRAMEWORK_DIR}" yes || failures=$((failures + 1))

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

compare_arms "${CASE_DIR}" "${FRAMEWORK_DIR}" "${FRAMEWORK_D_REL_TOL}" \
    "two materials, framework vs legacy" || failures=$((failures + 1))

compare_arms "${SINGLE_DIR}" "${SINGLE_FRAMEWORK_DIR}" "${SINGLE_D_REL_TOL}" \
    "one material, framework vs legacy" || failures=$((failures + 1))

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
    for dir in "${CASE_DIR}" "${FRAMEWORK_DIR}" "${SINGLE_DIR}" \
        "${SINGLE_FRAMEWORK_DIR}"
    do
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
