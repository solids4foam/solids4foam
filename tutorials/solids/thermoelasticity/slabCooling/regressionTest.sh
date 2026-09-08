#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# slabCooling regression test
# Unconstrained thermal contraction
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

# Stress should be ~0 (numerical noise only)
SIGMA_MAX=1e3      # Pa

# Strain should be O(1e-8)
EPS_MIN=1e-9
EPS_MAX=1e-7

# ------------------------------------------------------------
# Log files
# ------------------------------------------------------------

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

variant="openfoamcom"
if [[ -n "${FOAMEXTEND:-}" || "${WM_PROJECT_VERSION:-}" == "4.1" ]]; then
    variant="foamextend"
elif [[ "${WM_PROJECT_VERSION:-}" != *"v"* ]]; then
    variant="openfoamorg"
fi

if [[ "${variant}" != "openfoamcom" ]]; then
    SIGMA_MAX=1.05e3
fi

echo "============================================================"
echo "slabCooling regression test"
echo "Max sigmaEq < ${SIGMA_MAX} Pa"
echo "epsilonEq order: ${EPS_MIN} < eps < ${EPS_MAX}"
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
# Clean & run
# ------------------------------------------------------------

prepare_case
( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )

# Run the case a second time with the stress taken from the
# mechanicalConstitutiveLaw framework rather than the legacy mechanicalModel,
# and require the two to agree.
#
# thermoMechanicalLaw is the first composite law on the framework: it owns a
# sub-law, delegates to it, and subtracts the thermal term. The two arms differ
# in one dictionary entry and nothing else
run_framework_comparison() {
    local dir="${REGRESSION_ROOT}/framework"

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

    # The switch goes inside the solid model's coeffs block, which is where the
    # model looks for it; at the top level it is read by nothing
    sed -i \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${dir}/constant/solidProperties"

    ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
        echo "FAIL: the framework arm did not run"
        return 1
    }

    # Each arm must have taken the path it was set up for, or this compares two
    # copies of the same thing
    if ! grep -q "Selecting mechanical constitutive law" \
        "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the framework arm did not use the framework"
        return 1
    fi

    if grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the legacy arm used the framework"
        return 1
    fi

    local t
    t=$(foamListTimes -case "${dir}" -latestTime 2>/dev/null | tail -n 1)

    if [[ -z "${t}" || ! -f "${dir}/${t}/D" || ! -f "${CASE_DIR}/${t}/D" ]]; then
        echo "FAIL: the framework comparison produced no D field"
        return 1
    fi

    if diff -q "${CASE_DIR}/${t}/D" "${dir}/${t}/D" > /dev/null; then
        echo "PASS: framework and legacy agree exactly"
        return 0
    fi

    echo "FAIL: framework and legacy differ"
    return 1
}

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk -F '=' '{print $2}' \
        | tr -d '[:space:]'
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk -F '=' '{print $2}' \
        | tr -d '[:space:]'
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

epsilon=$(extract_max_epsilon)
sigma=$(extract_max_sigma)

if [[ -z "${epsilon}" || -z "${sigma}" ]]
then
    echo "FAIL: Could not extract epsilonEq or sigmaEq"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

# --- Stress check ------------------------------------------------------------

if awk "BEGIN {exit !(${sigma} < ${SIGMA_MAX})}"
then
    printf "PASS: Max sigmaEq = %.6g Pa\n" "${sigma}"
else
    printf "FAIL: Max sigmaEq = %.6g Pa\n" "${sigma}"
    failures=$((failures + 1))
fi

# --- Strain order-of-magnitude check ----------------------------------------

if awk "BEGIN {exit !(${epsilon} > ${EPS_MIN} && ${epsilon} < ${EPS_MAX})}"
then
    printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
else
    printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
    failures=$((failures + 1))
fi

# Clean case again
if ! run_framework_comparison; then
    failures=$((failures + 1))
fi

( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

echo
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
