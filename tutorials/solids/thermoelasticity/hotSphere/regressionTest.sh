#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

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

# Run the case a second time with the stress taken from the
# mechanicalConstitutiveLaw framework rather than the legacy mechanicalModel.
#
# This case runs several time steps, so unlike slabCooling it exercises the
# framework rolling its constitutive state over between them.
#
# The two arms are compared to a tolerance rather than exactly. They solve the
# same problem but reach it by slightly different iteration paths, because the
# implicit stiffness that steers the iteration is built differently: the
# framework interpolates its cell tangent to the faces where the legacy model
# forms a face value directly. The converged answers therefore agree only to
# the solution tolerance, and the difference shrinks with it - at
# solutionTolerance 1e-6 it is around 1e-7 of the displacement, and tightening
# to 1e-10 takes it to 1e-8. The threshold here is well above that and far
# below anything physical
FRAMEWORK_D_REL_TOL=1e-6

run_framework_comparison() {
    local legacy_dir="${REGRESSION_ROOT}/frameworkLegacy"
    local framework_dir="${REGRESSION_ROOT}/framework"
    local dir

    # Both arms are prepared and run here rather than reusing the main case,
    # so the comparison does not depend on what the main run left behind
    for dir in "${legacy_dir}" "${framework_dir}"; do
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
    done

    # Write enough digits for the comparison to be about the solution rather
    # than about the file format. At the default six significant figures the
    # two arms differ by around 3e-6 simply because that is the last digit
    # written, which would tell us nothing
    for dir in "${legacy_dir}" "${framework_dir}"; do
        if grep -q "^writePrecision" "${dir}/system/controlDict"; then
            sed -i 's|^writePrecision.*|writePrecision  14;|' \
                "${dir}/system/controlDict"
        else
            echo "writePrecision  14;" >> "${dir}/system/controlDict"
        fi
    done

    # The two arms differ in this one entry and nothing else. It goes inside
    # the solid model's coeffs block, which is where the model looks for it
    sed -i \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${framework_dir}/constant/solidProperties"

    for dir in "${legacy_dir}" "${framework_dir}"; do
        ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
            echo "FAIL: the framework comparison could not run ${dir}"
            return 1
        }
    done

    # Each arm must have taken the path it was set up for
    if ! grep -q "Selecting mechanical constitutive law" \
        "${framework_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the framework arm did not use the framework"
        return 1
    fi

    if grep -q "Selecting mechanical constitutive law" \
        "${legacy_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the legacy arm used the framework"
        return 1
    fi

    local tL tF
    tL=$(foamListTimes -case "${legacy_dir}" -latestTime 2>/dev/null | tail -n 1)
    tF=$(foamListTimes -case "${framework_dir}" -latestTime 2>/dev/null \
        | tail -n 1)

    if [[ -z "${tL}" || "${tL}" != "${tF}" ]]; then
        echo "FAIL: the two arms reached different times ('${tL}' vs '${tF}')"
        return 1
    fi

    if [[ ! -f "${legacy_dir}/${tL}/D" || ! -f "${framework_dir}/${tF}/D" ]]
    then
        echo "FAIL: the framework comparison produced no D field"
        return 1
    fi

    # Largest component-wise difference, relative to the largest component
    local rel
    rel=$(awk '
        function abs(x) { return x < 0 ? -x : x }
        FNR == NR {
            if ($0 ~ /^\(/)
            {
                gsub(/[()]/, ""); n++
                for (i = 1; i <= NF; i++) a[n, i] = $i
            }
            next
        }
        {
            if ($0 ~ /^\(/)
            {
                gsub(/[()]/, ""); k++
                for (i = 1; i <= NF; i++)
                {
                    d = abs($i - a[k, i]); if (d > maxd) maxd = d
                    v = abs(a[k, i]);      if (v > maxv) maxv = v
                }
            }
        }
        END { printf "%.10g\n", (maxv > 0 ? maxd/maxv : maxd) }
    ' "${legacy_dir}/${tL}/D" "${framework_dir}/${tF}/D")

    if [[ -z "${rel}" ]]; then
        echo "FAIL: could not compare the two arms"
        return 1
    fi

    if awk "BEGIN {exit !(${rel} < ${FRAMEWORK_D_REL_TOL})}"; then
        printf "PASS: framework and legacy agree, relative D diff = %.4g\n" \
            "${rel}"
        return 0
    fi

    printf "FAIL: framework and legacy differ, relative D diff = %.4g\n" \
        "${rel}"
    return 1
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
if ! run_framework_comparison; then
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
