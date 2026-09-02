#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
# Run twice: the legacy mechanicalModel, and the mechanicalConstitutiveLaw
# framework. This is the only regression case whose material is history
# dependent AND which runs the framework end to end in a solver, so it is the
# only place a plastic history error in the framework would show up in a real
# solve rather than in a unit check
APPROACHES=(
    legacy
    framework
)

# ============================================================
# Elastoplastic perforated plate regression test
# Checks strain, stress, and plastic yielding
# ============================================================

# Reference ranges (order-of-magnitude + robustness)
EPSILON_MIN=4.5e-3
EPSILON_MAX=6.0e-3

SIGMA_MIN=1.5e8
SIGMA_MAX=2.2e8

YIELD_MIN=28
YIELD_MAX=32

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Elastoplastic perforated plate regression test"
echo "Max epsilonEq           in [${EPSILON_MIN}, ${EPSILON_MAX}]"
echo "Max sigmaEq (von Mises) in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "Yielding cells          in [${YIELD_MIN}, ${YIELD_MAX}]"
echo "============================================================"
echo

prepare_case() {
    local approach="$1"
    CASE_DIR="${REGRESSION_ROOT}/${approach}"

    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done

    if [[ "${approach}" == "framework" ]]; then
        # The switch lives in the <type>Coeffs sub-dictionary; at the top level
        # it is silently ignored and this arm would quietly repeat the legacy run
        sed -i.bak \
            's/^    predictor yes;/    useMechanicalConstitutiveLawManager yes;\n    predictor yes;/' \
            "${CASE_DIR}/constant/solidProperties"
        rm -f "${CASE_DIR}/constant/solidProperties.bak"

        if ! grep -q 'useMechanicalConstitutiveLawManager' \
            "${CASE_DIR}/constant/solidProperties"
        then
            echo "FAIL: could not enable the framework in solidProperties"
            exit 1
        fi
    fi
}

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

extract_yielding_cells() {
    # The legacy law counts cells; the framework law counts integration points,
    # which is the honest description for a law that may be evaluated on faces
    # or points. Accept either wording rather than making one law lie
    # Prefer the framework's message: in the framework arm the legacy law is
    # still constructed but never called, so it truthfully reports zero and
    # would mask the framework's own count
    if grep -q "yielding integration points" "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        grep "Number of yielding integration points" \
            "${CASE_DIR}/${SOLVER_LOGFILE}" \
            | tail -n 101 \
            | head -n 1 \
            | sed 's|.*= *||; s|/.*||'
    elif grep -q "cells .* are actively yielding" "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        grep "cells .* are actively yielding" "${CASE_DIR}/${SOLVER_LOGFILE}" \
            | tail -n 101 \
            | head -n 1 \
            | awk '{print $1}'
    fi
}

# Exercise the mechanicalConstitutiveLawManager on this case. It is the only
# tutorial in the regression set whose material is history dependent, so it is
# the only one where a tangent query has any history to disturb
run_constitutive_test() {
    if ! command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
        echo "SKIP: Test-mechanicalConstitutiveLaw not found in PATH"
        return 0
    fi

    if [[ ! -d "${CASE_DIR}/constant/polyMesh" ]]; then
        echo "SKIP: mechanicalConstitutiveLaw checks (case has no mesh)"
        return 0
    fi

    if ( cd "${CASE_DIR}" && Test-mechanicalConstitutiveLaw \
            > "log.Test-mechanicalConstitutiveLaw" 2>&1 )
    then
        local n_passed
        n_passed=$(grep -c 'PASS:' \
            "${CASE_DIR}/log.Test-mechanicalConstitutiveLaw" || true)

        if (( n_passed == 0 )); then
            echo "SKIP: mechanicalConstitutiveLaw checks (not built on this fork)"
            return 0
        fi

        echo "PASS: mechanicalConstitutiveLaw checks (${n_passed} checks)"
        return 0
    fi

    echo "FAIL: mechanicalConstitutiveLaw checks"
    grep 'FAIL:' "${CASE_DIR}/log.Test-mechanicalConstitutiveLaw" || true
    return 1
}

# ------------------------------------------------------------
# Run and check each approach
# ------------------------------------------------------------

failures=0
declare -A RESULT_EPS
declare -A RESULT_SIG

for approach in "${APPROACHES[@]}"; do
    echo
    echo "------------------------------------------------------------"
    echo "Testing approach: ${approach}"
    echo "------------------------------------------------------------"

    prepare_case "${approach}"

    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )

    # Assert the run took the path this approach names
    if grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        used_framework=true
    else
        used_framework=false
    fi

    if [[ "${approach}" == "framework" && "${used_framework}" == false ]]; then
        echo "FAIL: framework approach did not construct the framework"
        failures=$((failures + 1))
    elif [[ "${approach}" == "legacy" && "${used_framework}" == true ]]; then
        echo "FAIL: legacy approach unexpectedly constructed the framework"
        failures=$((failures + 1))
    else
        echo "PASS: ${approach} took the expected path"
    fi

    epsilon=$(extract_max_epsilon)
    sigma=$(extract_max_sigma)
    yielding_cells=$(extract_yielding_cells)

    if [[ -z "${epsilon}" || -z "${sigma}" || -z "${yielding_cells}" ]]; then
        echo "FAIL: ${approach} could not extract one or more quantities"
        failures=$((failures + 1))
        continue
    fi

    RESULT_EPS["${approach}"]="${epsilon}"
    RESULT_SIG["${approach}"]="${sigma}"

    if awk "BEGIN {exit !(${epsilon} >= ${EPSILON_MIN} && ${epsilon} <= ${EPSILON_MAX})}"
    then
        printf "PASS: %s Max epsilonEq = %.6g\n" "${approach}" "${epsilon}"
    else
        printf "FAIL: %s Max epsilonEq = %.6g\n" "${approach}" "${epsilon}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${sigma} <= ${SIGMA_MAX})}"
    then
        printf "PASS: %s Max sigmaEq = %.6g\n" "${approach}" "${sigma}"
    else
        printf "FAIL: %s Max sigmaEq = %.6g\n" "${approach}" "${sigma}"
        failures=$((failures + 1))
    fi

    if (( yielding_cells >= YIELD_MIN && yielding_cells <= YIELD_MAX )); then
        printf "PASS: %s yielding = %d\n" "${approach}" "${yielding_cells}"
    else
        printf "FAIL: %s yielding = %d\n" "${approach}" "${yielding_cells}"
        failures=$((failures + 1))
    fi

    if [[ "${approach}" == "legacy" ]]; then
        if ! run_constitutive_test; then
            failures=$((failures + 1))
        fi
    fi

    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
done

# The framework must reproduce the legacy result. This case has no pressure
# smoothing, so unlike longWall there is nothing the framework omits: the two
# should agree to solver tolerance, and a plastic history error would show
# here as a difference that the unit checks cannot see
if [[ -n "${RESULT_EPS[legacy]:-}" && -n "${RESULT_EPS[framework]:-}" ]]; then
    for q in eps sig; do
        if [[ "${q}" == "eps" ]]; then
            a="${RESULT_EPS[legacy]}"; b="${RESULT_EPS[framework]}"; n="epsilonEq"
        else
            a="${RESULT_SIG[legacy]}"; b="${RESULT_SIG[framework]}"; n="sigmaEq"
        fi

        if awk "BEGIN {exit !(($a - $b)^2 <= (1e-8*$a)^2)}"; then
            printf "PASS: legacy and framework %s agree (%.8g vs %.8g)\n" \
                "$n" "$a" "$b"
        else
            printf "FAIL: legacy and framework %s differ (%.8g vs %.8g)\n" \
                "$n" "$a" "$b"
            failures=$((failures + 1))
        fi
    done
else
    echo "SKIP: cross-check needs both approaches to have run"
fi

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
