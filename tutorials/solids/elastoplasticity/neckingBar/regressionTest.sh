#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
# Run twice: the legacy mechanicalModel and the mechanicalConstitutiveLaw
# framework. This is the end-to-end check of two things nothing else covers:
# the stress migration of nonLinGeomUpdatedLagSolid, and the finite-strain
# plastic law. It is the only tutorial pairing them that does not enable
# pressure smoothing, so it is the only one where the framework is expected to
# reproduce the legacy result exactly rather than approximately
APPROACHES=(
    legacy
    framework
)

# ============================================================
# neckingBar regression test
# Checks the final loading force after the necking curve.
# ============================================================

FORCE_MIN=74.3
FORCE_MAX=74.6

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "neckingBar regression test"
echo "Final loading force in [${FORCE_MIN}, ${FORCE_MAX}]"
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
        # it is silently ignored and this arm would repeat the legacy run
        sed -i.bak \
            's/^    nCorrectors     1000;/    useMechanicalConstitutiveLawManager yes;\n    nCorrectors     1000;/' \
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

failures=0
declare -A RESULT_F

for approach in "${APPROACHES[@]}"; do
    echo
    echo "------------------------------------------------------------"
    echo "Testing approach: ${approach}"
    echo "------------------------------------------------------------"

    if [ "$CHECK_ONLY" = false ]; then
        prepare_case "${approach}"
        sed -i.bak 's/^endTime         1;/endTime         0.184;/' \
            "${CASE_DIR}/system/controlDict"
        rm -f "${CASE_DIR}/system/controlDict.bak"
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
        ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
    else
        CASE_DIR="${REGRESSION_ROOT}/${approach}"
        echo "Running in check-only mode: skipping Allclean and Allrun"
    fi

    if grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/log.solids4Foam" 2>/dev/null
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

    force_file=$(find "${CASE_DIR}/postProcessing" \
        -name 'solidForcesDisplacementsloading.dat' -print | tail -n 1)
    if [[ -z "${force_file}" ]]; then
        echo "FAIL: ${approach} could not find loading force history"
        failures=$((failures + 1))
        continue
    fi

    final_force=$(awk 'END {print $5*360e-3}' "${force_file}")
    RESULT_F["${approach}"]="${final_force}"

    if awk "BEGIN {exit !(${final_force} >= ${FORCE_MIN} && ${final_force} <= ${FORCE_MAX})}"; then
        printf "PASS: %s final loading force = %.6g\n" "${approach}" "${final_force}"
    else
        printf "FAIL: %s final loading force = %.6g\n" "${approach}" "${final_force}"
        failures=$((failures + 1))
    fi

    if [ "$CHECK_ONLY" = false ]; then
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    fi
done

# No pressure smoothing here, so the framework omits nothing and the two must
# agree. A history error in the finite-strain plastic path, or in the relative
# deformation gradient the framework is given, shows up here and nowhere else
if [[ -n "${RESULT_F[legacy]:-}" && -n "${RESULT_F[framework]:-}" ]]; then
    a="${RESULT_F[legacy]}"; b="${RESULT_F[framework]}"

    if awk "BEGIN {exit !(($a - $b)^2 <= (1e-6*$a)^2)}"; then
        printf "PASS: legacy and framework force agree (%.8g vs %.8g)\n" "$a" "$b"
    else
        printf "FAIL: legacy and framework force differ (%.8g vs %.8g)\n" "$a" "$b"
        failures=$((failures + 1))
    fi
else
    echo "SKIP: cross-check needs both approaches to have run"
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
