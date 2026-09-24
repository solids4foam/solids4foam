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

# The final D of the removed legacy mechanicalModel, as a digest of its
# internal values as written, from the last commit that had it
# (mcl-stage8-coverage, c3a92b3d), per fork. thermoMechanicalLaw is a
# composite: it owns a sub-law, delegates to it, and subtracts the thermal
# term, and the framework reproduced the legacy D field exactly, in every
# figure written
case "$(solids4Foam::foamFlavour)" in
    com)
        LEGACY_D_DIGEST=b7609301dbc9037471029473cdc64be52c014832
        ;;
    org)
        LEGACY_D_DIGEST=b7609301dbc9037471029473cdc64be52c014832
        ;;
    foamextend)
        LEGACY_D_DIGEST=4760637f46f33deb124c59b8219f9e0a0d4e014e
        ;;
esac

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

# A digest of a field's internal values as written, independent of the file
# header, which names the OpenFOAM version that wrote it
internal_field_digest() {
    python3 - "$1" << 'PYEOF'
import hashlib
import re
import sys

text = open(sys.argv[1]).read()
field = re.search(r"\binternalField\s+(.*?);\s*\n\s*boundaryField", text, re.DOTALL)
if not field:
    sys.exit(f"cannot find the internalField in {sys.argv[1]}")
print(hashlib.sha1(" ".join(field.group(1).split()).encode()).hexdigest())
PYEOF
}

check_against_legacy() {
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

    local digest
    digest=$(internal_field_digest "${CASE_DIR}/${t}/D" || true)

    if [[ -n "${digest}" && "${digest}" == "${LEGACY_D_DIGEST}" ]]; then
        echo "PASS: D matches the legacy model's exactly"
        return 0
    fi

    echo "FAIL: D differs from the legacy model's" \
        "(digest ${digest:-none}, legacy ${LEGACY_D_DIGEST})"
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

if ! check_against_legacy; then
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
