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
# immersedHronTurekFsi2 regression test
# Runs the coarsest mesh (MESH_LEVEL=1) with the coupling started at
# t = 0.5 s rather than 2 s, to t = 0.6 s (100 coupled time steps), and
# checks the tip displacement of the flag and the force on it from the
# surface traction at the end time against the values of this method on
# this mesh.
# ============================================================

REG_COUPLING_START_TIME=0.5
REG_END_TIME=0.6

# Reference values at REG_END_TIME and tolerances
REF_TIP_UX=-2.00967e-05
REF_TIP_UY=-0.000138216
REF_FX=-0.258839
REF_FY=0.0691452
UX_TOL=2e-6
UY_TOL=5e-6
FX_TOL=5e-3
FY_TOL=5e-3

ALLRUN_LOGFILE="log.Allrun"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"
FORCE_FILE="postProcessing/fluid/immersedBoundary/0/flag.dat"

echo "============================================================"
echo "immersedHronTurekFsi2 regression test"
echo "Coupling start time = ${REG_COUPLING_START_TIME}"
echo "Regression end time = ${REG_END_TIME}"
echo "Tip Ux tolerance    < ${UX_TOL}"
echo "Tip Uy tolerance    < ${UY_TOL}"
echo "Fx tolerance        < ${FX_TOL}"
echo "Fy tolerance        < ${FY_TOL}"
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
    sed -i.bak \
        "s/^\([[:space:]]*couplingStartTime[[:space:]]\+\)[0-9.]*;/\1${REG_COUPLING_START_TIME};/" \
        "${CASE_DIR}/constant/fsiProperties"
    rm -f "${CASE_DIR}/system/controlDict.bak" \
        "${CASE_DIR}/constant/fsiProperties.bak"
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

for f in "${DISP_FILE}" "${FORCE_FILE}"; do
    if [[ ! -f "${CASE_DIR}/${f}" ]]; then
        echo "FAIL: Could not find ${f}"
        exit 1
    fi
done

# Tip displacement at the end time
IFS=" " read -r tip_time tip_ux tip_uy < <(awk '
    ($1 + 0) == $1 { t = $1; ux = $2; uy = $3 }
    END { if (t != "") print t, ux, uy }
' "${CASE_DIR}/${DISP_FILE}")

# Force on the flag from the surface traction (columns 11-13) at the end time
IFS=" " read -r force_time fx fy < <(awk '
    !/^#/ { t = $1; fx = $11; fy = $12 }
    END { if (t != "") print t, fx, fy }
' "${CASE_DIR}/${FORCE_FILE}")

if [[ -z "${tip_time:-}" || -z "${force_time:-}" ]]; then
    echo "FAIL: Could not extract the tip displacement or the force"
    exit 1
fi

if ! awk "BEGIN {exit !(${tip_time} + 0 >= ${REG_END_TIME} - 1e-6)}"; then
    echo "FAIL: the tip displacement history stops at t = ${tip_time}, short of"
    echo "      the end time ${REG_END_TIME}: the case did not complete"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check_value() {
    local name=$1 value=$2 ref=$3 tol=$4
    if awk "BEGIN {d = ${value} - ${ref}; if (d < 0) d = -d; exit !(d < ${tol})}"; then
        printf "PASS: %s = %.6g (reference %.6g)\n" "${name}" "${value}" "${ref}"
    else
        printf "FAIL: %s = %.6g (reference %.6g, tolerance %g)\n" \
            "${name}" "${value}" "${ref}" "${tol}"
        failures=$((failures + 1))
    fi
}

check_value "Tip Ux" "${tip_ux}" "${REF_TIP_UX}" "${UX_TOL}"
check_value "Tip Uy" "${tip_uy}" "${REF_TIP_UY}" "${UY_TOL}"
check_value "Fx" "${fx}" "${REF_FX}" "${FX_TOL}"
check_value "Fy" "${fy}" "${REF_FY}" "${FY_TOL}"

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
