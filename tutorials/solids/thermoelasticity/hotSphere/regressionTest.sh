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
# hotSphere regression test
# Uses the tutorial's reported temperature and stress extrema.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

T_MIN=339.5
T_MAX=340.5
EPS_MIN=1.5e-4
EPS_MAX=2.5e-4
SIGMA_MIN=3.8e7
SIGMA_MAX=4.6e7

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "hotSphere regression test"
echo "Max T       in [${T_MIN}, ${T_MAX}]"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq in [${SIGMA_MIN}, ${SIGMA_MAX}]"
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
}

shorten_case() {
    local controlDict="${CASE_DIR}/system/controlDict"
    sed -i.bak 's/^endTime[[:space:]]\+5;/endTime         1;/' "${controlDict}"
    rm -f "${controlDict}.bak"
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
    shorten_case
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_temperature() {
    grep "Max T" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk '{print $NF}'
}

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk '{print $NF}'
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk '{print $NF}'
}

# Run the case to its full end time, and hold it to the answer of the removed
# legacy mechanicalModel.
#
# This case runs several time steps, so unlike slabCooling it exercises the
# framework rolling its constitutive state over between them.
#
# The comparison is to a tolerance rather than exact. The framework and the
# legacy model solved the same problem but reached it by slightly different
# iteration paths, because the implicit stiffness that steers the iteration is
# built differently: the framework interpolates its cell tangent to the faces
# where the legacy model formed a face value directly. The converged answers
# therefore agreed only to the solution tolerance, and the difference shrank
# with it - at solutionTolerance 1e-6 it was around 1e-7 of the displacement,
# and tightening to 1e-10 took it to 1e-8. The threshold here is well above
# that and far below anything physical
FRAMEWORK_D_REL_TOL=1e-6
COMPARISON_END_TIME=5

# The legacy model's final D, as the max and mean component magnitude of the
# field written to fourteen figures, from the last commit that had it
# (mcl-stage8-coverage, c3a92b3d), per fork
case "$(solids4Foam::foamFlavour)" in
    com)
        LEGACY_D_MAX=0.00015702498751748
        LEGACY_D_MEAN=7.86486280639836e-05
        ;;
    org)
        LEGACY_D_MAX=0.00015702498782177
        LEGACY_D_MEAN=7.86486272668098e-05
        ;;
    *)
        # The case does not run here
        LEGACY_D_MAX=""
        LEGACY_D_MEAN=""
        ;;
esac

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

run_legacy_comparison() {
    # Not on foam-extend, where the framework and the legacy model differed by
    # 0.6 % in D, and the framework's answer is the one kept. They matched to
    # 1e-13 for two correctors and part of the third, and differed only in the
    # stress on the three symmetryPlane patches. The framework corrects the
    # stress's boundary conditions after evaluating it; the legacy law did not.
    # On foam-extend's symmetryPlane that correction replaces the law's
    # boundary value with the symmetry transform of the adjacent cell's
    # stress, which has zero shear traction on the plane, as symmetry
    # requires. Correcting the legacy stress too reproduced the framework's
    # answer to 2e-8.
    #
    # The transform was checked because older OpenFOAM versions were thought to
    # get it wrong for tensors: 0.5*(x + transform(I - 2nn, x)) matches a
    # hand-written reference to 4e-16 for symmTensor and tensor, with zero shear
    # traction and all other components kept, on v1912-v2606, OpenFOAM.org 8-13
    # and foam-extend 4.1 and 5.0. So the legacy answer there is the one known
    # to be wrong, and there is nothing to hold the framework to
    if [[ "${WM_PROJECT:-}" == "foam" ]]; then
        echo "SKIP: legacy comparison (known foam-extend symmetryPlane difference)"
        return 0
    fi

    local dir="${REGRESSION_ROOT}/comparison"

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

    # Write enough digits for the comparison to be about the solution rather
    # than about the file format. At the default six significant figures the
    # two differ by around 3e-6 simply because that is the last digit written,
    # which would tell us nothing
    if grep -q "^writePrecision" "${dir}/system/controlDict"; then
        sed -i 's|^writePrecision.*|writePrecision  14;|' \
            "${dir}/system/controlDict"
    else
        echo "writePrecision  14;" >> "${dir}/system/controlDict"
    fi

    ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
        echo "FAIL: the legacy comparison could not run ${dir}"
        return 1
    }

    if ! grep -q "Selecting mechanical constitutive law" \
        "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the comparison constructed no mechanical constitutive law"
        return 1
    fi

    if ! grep -q "^End" "${dir}/${SOLVER_LOGFILE}" \
      || grep -qE "Nonlinear solve did not converge|SNES convergence error|FOAM FATAL" \
          "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${dir} did not complete and converge"
        return 1
    fi

    local t
    t=$(solids4Foam::latestTime "${dir}")

    if [[ -z "${t}" ]] \
        || ! awk "BEGIN {exit !((${t} - ${COMPARISON_END_TIME})^2 <= 1e-20)}"
    then
        echo "FAIL: comparison stopped at '${t}'; expected ${COMPARISON_END_TIME}"
        return 1
    fi

    check_field_against_legacy "D at t = ${t}" "${dir}/${t}/D" \
        "${LEGACY_D_MAX}" "${LEGACY_D_MEAN}" "${FRAMEWORK_D_REL_TOL}"
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

temperature=$(extract_max_temperature)
epsilon=$(extract_max_epsilon)
sigma=$(extract_max_sigma)

if [[ -z "${temperature}" || -z "${epsilon}" || -z "${sigma}" ]]; then
    echo "FAIL: Could not extract one or more regression quantities"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${temperature} >= ${T_MIN} && ${temperature} <= ${T_MAX})}"; then
    printf "PASS: Max T = %.6g\n" "${temperature}"
else
    printf "FAIL: Max T = %.6g\n" "${temperature}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${epsilon} >= ${EPS_MIN} && ${epsilon} <= ${EPS_MAX})}"; then
    printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
else
    printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${SIGMA_MAX} >= ${sigma})}"; then
    printf "PASS: Max sigmaEq = %.6g\n" "${sigma}"
else
    printf "FAIL: Max sigmaEq = %.6g\n" "${sigma}"
    failures=$((failures + 1))
fi

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
    if ! run_legacy_comparison; then
        failures=$((failures + 1))
    fi

    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
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
