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
# Rod and seabed regression test
# Checks strain and stress
# ============================================================

# Reference ranges (order-of-magnitude + robustness)
#
# These moved when the anisotropic Biot law stopped selecting its reduced
# plane model on this three-dimensional mesh. The old values were not
# self-consistent with the declared material: 80 kPa against a strain of
# 9.3e-5 implies a stiffness near 9e8 Pa, where the moduli here are 1.2e7 to
# 2e7 Pa. Forcing the out-of-plane stress to zero while xx and yy were not
# manufactured a large deviator, so von Mises read high against a small
# strain. The values below sit on the declared moduli
EPSILON_MIN=1.4e-3
EPSILON_MAX=2.2e-3

SIGMA_MIN=40e3
SIGMA_MAX=58e3

# The final D of the removed legacy mechanicalModel, as the max and mean
# component magnitude of the field written to fourteen figures, from the last
# commit that had it (mcl-stage8-coverage, c3a92b3d), per fork. The case as it
# ships is poroMechanicalLaw over anisotropicBiotElastic, and the framework
# reproduced the legacy D field exactly, in every one of those figures. These are
# recorded numbers, though, and another compiler, CPU or MPI build moves an
# iterative solution by round-off at the solver tolerance: CI measures up to
# 3e-8 relative against values recorded on macOS. The tolerance, 1e-6 of the
# largest value, allows for that, and is ten times below the 1e-5 that a
# 0.001% change in a material constant makes
case "$(solids4Foam::foamFlavour)" in
    com)
        LEGACY_D_MAX=0.024504218766989
        LEGACY_D_MEAN=0.00528889093233214
        ;;
    org)
        LEGACY_D_MAX=0.024504219027988
        LEGACY_D_MEAN=0.00528889109572662
        ;;
    foamextend)
        LEGACY_D_MAX=0.024496819721199
        LEGACY_D_MEAN=0.00528997493008025
        ;;
esac
LEGACY_D_REL_TOL=1e-6

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Road and seabed regression test"
echo "Max epsilonEq           in [${EPSILON_MIN}, ${EPSILON_MAX}]"
echo "Max sigmaEq (von Mises) in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done

    # Enough digits that the comparison with the legacy answer is about the
    # solution and not about the last figure written
    if grep -q "^writePrecision" "${CASE_DIR}/system/controlDict"; then
        sed -i 's|^writePrecision.*|writePrecision  14;|' \
            "${CASE_DIR}/system/controlDict"
    else
        echo "writePrecision  14;" >> "${CASE_DIR}/system/controlDict"
    fi
}

# ------------------------------------------------------------
# Clean & run case
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

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

epsilon=$(extract_max_epsilon)
sigma=$(extract_max_sigma)

if [[ -z "${epsilon}" || -z "${sigma}" ]]
then
    echo "FAIL: Could not extract one or more regression quantities"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

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

# Check the poroMechanicalLaw composite against the legacy law, on the case as
# it ships: poroMechanicalLaw over anisotropicBiotElastic.
#
# This is the case that exercises the effective stress the composite carries.
# anisotropicBiotElastic leaves the zz, yz and xz components of the stress
# unwritten in the branch this case takes, so they come from whatever the
# sub-law was given to work in - which is the whole reason the composite hands
# it the effective stress rather than the caller's total stress
check_poro_against_legacy() {
    if ! grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the case constructed no mechanical constitutive law"
        return 1
    fi

    local t end_time
    t=$(solids4Foam::latestTime "${CASE_DIR}")
    end_time=$(sed -n 's/^endTime[[:space:]]*\([^;]*\);.*/\1/p' \
        "${CASE_DIR}/system/controlDict")

    if [[ -z "${t}" || -z "${end_time}" ]] \
        || ! awk "BEGIN {exit !((${t} - ${end_time})^2 <= 1e-20)}"
    then
        echo "FAIL: the case stopped at '${t}', not at the end time '${end_time}'"
        return 1
    fi

    check_field_against_legacy "poro D at t = ${t}" "${CASE_DIR}/${t}/D" \
        "${LEGACY_D_MAX}" "${LEGACY_D_MEAN}" "${LEGACY_D_REL_TOL}"
}

failures=0

# --- epsilonEq ---
if awk "BEGIN {exit !(${epsilon} >= ${EPSILON_MIN} && ${epsilon} <= ${EPSILON_MAX})}"
then
    printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
else
    printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
    failures=$((failures + 1))
fi

# --- sigmaEq ---
if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${sigma} <= ${SIGMA_MAX})}"
then
    printf "PASS: Max sigmaEq = %.6g\n" "${sigma}"
else
    printf "FAIL: Max sigmaEq = %.6g\n" "${sigma}"
    failures=$((failures + 1))
fi

echo
if ! check_poro_against_legacy; then
    failures=$((failures + 1))
fi

if (( failures == 0 ))
then
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
