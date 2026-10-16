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

# GNU sed, for the in-place edits below
solids4Foam::requireGnuSed

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

# The final D, as the max and mean component magnitude of the field written
# to fourteen figures. The case as it ships is poroMechanicalLaw over
# anisotropicBiotElastic, and the framework reproduced the removed legacy
# mechanicalModel's D field exactly, in every one of those figures, on every
# fork (mcl-stage8-coverage, c3a92b3d). These are OpenFOAM.com v2512's;
# OpenFOAM.org 9 agrees to 1e-8, and foam-extend 4.1 gives a max 3e-4 smaller.
# The tolerance, 5e-4 of the largest value, covers that
REF_D_MAX=0.024504218766989
REF_D_MEAN=0.00528889093233214
REF_D_REL_TOL=5e-4

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
        "${SOLIDS4FOAM_SED}" -i 's|^writePrecision.*|writePrecision  14;|' \
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

# Check the poroMechanicalLaw composite against the legacy law, on the case as
# it ships: poroMechanicalLaw over anisotropicBiotElastic.
#
# This is the case that exercises the effective stress the composite carries.
# anisotropicBiotElastic leaves the zz, yz and xz components of the stress
# unwritten in the branch this case takes, so they come from whatever the
# sub-law was given to work in - which is the whole reason the composite hands
# it the effective stress rather than the caller's total stress
check_poro_against_reference() {
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

    solids4Foam::checkFieldNorms "poro D at t = ${t}" "${CASE_DIR}/${t}/D" \
        "${REF_D_MAX}" "${REF_D_MEAN}" "${REF_D_REL_TOL}"
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
if ! check_poro_against_reference; then
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
