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
# bendingBeamInChannel regression test
# Runs the coarsest mesh (MESH_LEVEL=1) to t = 2 s and checks
# the mean and root mean square drag coefficient of the
# bending immersed beam over its second cycle, 1 < t < 2 s,
# against the values of this method on this mesh.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

REG_END_TIME=2
CD_MEAN_MIN=17.7
CD_MEAN_MAX=18.7
CD_RMS_MIN=53.8
CD_RMS_MAX=56.0

ALLRUN_LOGFILE="log.Allrun"
FORCE_FILE="postProcessing/immersedBoundary/0/beam.dat"

echo "============================================================"
echo "bendingBeamInChannel regression test"
echo "Regression end time = ${REG_END_TIME}"
echo "Mean Cd in [${CD_MEAN_MIN}, ${CD_MEAN_MAX}]"
echo "RMS Cd in [${CD_RMS_MIN}, ${CD_RMS_MAX}]"
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

# Cd = 2 Fx/(rho Umean^2 H Lz) = 9.855 Fx, with the force from the surface
# traction (columns 2-4), time-averaged over 1 < t < 2 s
read -r mean_cd rms_cd < <(awk '
    !/^#/ {
        if (tOld != "" && $1 > 1) {
            dt = $1 - tOld; cd = 9.855*$2
            s += cd*dt; s2 += cd*cd*dt; T += dt
        }
        tOld = $1
    }
    END { if (T > 0) printf "%.6g\t%.6g\n", s/T, sqrt(s2/T) }
' "${CASE_DIR}/${FORCE_FILE}")

if [[ -z "${mean_cd:-}" || -z "${rms_cd:-}" ]]; then
    echo "FAIL: Could not extract the drag coefficients over 1 < t < 2 s"
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

check_range "Mean Cd" "${mean_cd}" "${CD_MEAN_MIN}" "${CD_MEAN_MAX}"
check_range "RMS Cd" "${rms_cd}" "${CD_RMS_MIN}" "${CD_RMS_MAX}"

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
