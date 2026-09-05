#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

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

# Check the poroMechanicalLaw composite against the legacy law, on the case as
# it ships: poroMechanicalLaw over anisotropicBiotElastic. The two arms differ
# in one dictionary entry and nothing else, and must agree exactly.
#
# This is the case that exercises the effective stress the composite carries.
# anisotropicBiotElastic leaves the zz, yz and xz components of the stress
# unwritten in the branch this case takes, so they come from whatever the
# sub-law was given to work in - which is the whole reason the composite hands
# it the effective stress rather than the caller's total stress
run_poro_framework_comparison() {
    local legacy_dir="${REGRESSION_ROOT}/poroLegacy"
    local framework_dir="${REGRESSION_ROOT}/poroFramework"
    local dir

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

        # Enough digits that the comparison is about the solution and not
        # about the last figure written
        if grep -q "^writePrecision" "${dir}/system/controlDict"; then
            sed -i 's|^writePrecision.*|writePrecision  14;|' \
                "${dir}/system/controlDict"
        else
            echo "writePrecision  14;" >> "${dir}/system/controlDict"
        fi
    done

    sed -i \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${framework_dir}/constant/solidProperties"

    for dir in "${legacy_dir}" "${framework_dir}"; do
        ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
            echo "FAIL: the poro comparison could not run ${dir}"
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
        echo "FAIL: the poro arms reached different times ('${tL}' vs '${tF}')"
        return 1
    fi

    if [[ ! -f "${legacy_dir}/${tL}/D" || ! -f "${framework_dir}/${tF}/D" ]]
    then
        echo "FAIL: the poro comparison produced no D field"
        return 1
    fi

    if diff -q "${legacy_dir}/${tL}/D" "${framework_dir}/${tF}/D" > /dev/null
    then
        echo "PASS: poro framework and legacy agree exactly"
        return 0
    fi

    echo "FAIL: poro framework and legacy differ"
    return 1
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
if ! run_poro_framework_comparison; then
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
