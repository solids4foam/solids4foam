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
# oscillatingCylinderInChannel regression test
# Runs the coarsest mesh (MESH_LEVEL=1) for a quarter of a
# period and checks the drag and lift coefficients of the
# immersed oscillating cylinder against the values of this
# method on this mesh.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

REG_END_TIME=1
CD_MIN=1.25
CD_MAX=1.35
CL_MIN=-0.015
CL_MAX=-0.002

# Opt-in check of the clustered surface-traction quadrature:
#     TRACTION_QUADRATURE=clustered ./regressionTest.sh
# The tolerances above are those of the default (greedy) quadrature, whose
# lift has a drag-to-lift leak that depends on the order of the STL facets.
# With the clustered quadrature the lift at t = 1 is that of the quadrature
# with every candidate point (0.0100 on this mesh), and that of the
# body-fitted reference solutions of verification/reference is 0.0095
# (CdBodyFitted/ClBodyFitted) to 0.0097 (bodyFitted901k.csv): the range
# covers them with a margin of about 20%. The drag is unchanged
TRACTION_QUADRATURE="${TRACTION_QUADRATURE:-greedy}"
if [[ "${TRACTION_QUADRATURE}" == "clustered" ]]; then
    CL_MIN=0.008
    CL_MAX=0.012
elif [[ "${TRACTION_QUADRATURE}" != "greedy" ]]; then
    echo "FAIL: TRACTION_QUADRATURE must be greedy (default) or clustered"
    exit 1
fi

ALLRUN_LOGFILE="log.Allrun"
FORCE_FILE="postProcessing/immersedBoundary/0/cylinder.dat"

echo "============================================================"
echo "oscillatingCylinderInChannel regression test"
echo "Regression end time = ${REG_END_TIME}"
if [[ "${TRACTION_QUADRATURE}" == "clustered" ]]; then
    echo "Traction quadrature = clustered (opt-in check)"
fi
echo "Final Cd in [${CD_MIN}, ${CD_MAX}]"
echo "Final Cl in [${CL_MIN}, ${CL_MAX}]"
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

    sed -i.bak \
        "s/^endTime[[:space:]]\+[0-9.]*;/endTime         ${REG_END_TIME};/" \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"

    if [[ "${TRACTION_QUADRATURE}" == "clustered" ]]; then
        local fvOptions="${CASE_DIR}/constant/fvOptions"
        sed -i.bak \
            "/type *immersedBoundaryForce;/a\\    tractionQuadrature clustered;" \
            "${fvOptions}"
        rm -f "${fvOptions}.bak"
        if ! grep -q "tractionQuadrature clustered;" "${fvOptions}"; then
            echo "FAIL: could not select tractionQuadrature clustered"
            exit 1
        fi
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
    ( cd "${CASE_DIR}" && MESH_LEVEL=1 ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
    echo "FAIL: Could not find ${FORCE_FILE}"
    exit 1
fi

# Cd = 2 Fx/(rho Uref^2 D Lz) = 1296.9 Fx, and likewise for Cl, from the force
# from the surface traction (forceEstimator surfaceTraction), which includes
# the inertia of the fluid inside the cylinder
final_cd=$(awk '!/^#/ { cd = 1296.9*$2 } END { print cd }' "${CASE_DIR}/${FORCE_FILE}")
final_cl=$(awk '!/^#/ { cl = 1296.9*$3 } END { print cl }' "${CASE_DIR}/${FORCE_FILE}")

if [[ -z "${final_cd}" || -z "${final_cl}" ]]; then
    echo "FAIL: Could not extract the final drag and lift coefficients"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check_range() {
    local name=$1 value=$2 min=$3 max=$4
    if awk "BEGIN {exit !(${value} >= ${min} && ${value} <= ${max})}"; then
        printf "PASS: %s = %.6g\n" "${name}" "${value}"
    else
        printf "FAIL: %s = %.6g\n" "${name}" "${value}"
        failures=$((failures + 1))
    fi
}

check_range "Final Cd" "${final_cd}" "${CD_MIN}" "${CD_MAX}"
check_range "Final Cl" "${final_cl}" "${CL_MIN}" "${CL_MAX}"

# Clean case again
if [ "$CHECK_ONLY" = false ]; then
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
