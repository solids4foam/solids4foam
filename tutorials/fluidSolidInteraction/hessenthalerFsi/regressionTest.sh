#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# hessenthalerFsi FSI regression test
#
# The full Phase I case runs for many seconds of simulated time on many
# cores, so the regression test runs only the first few time steps of the
# ramp-up, on the coarse fluid mesh, in a copy of the case.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances and reference values at REG_END_TIME
# ------------------------------------------------------------

# Five time steps of 2e-3 s
REG_END_TIME=0.01

# Flap tip displacement components (m)
REF_TIP_DY=-1.84381e-07
REF_TIP_DZ=-2.35088e-08
TIP_TOL=2e-08

ALLRUN_LOGFILE="log.Allrun"
TIP_FILE="postProcessing/0/solidPointDisplacement_flapTip.dat"

echo "============================================================"
echo "hessenthalerFsi FSI regression test"
echo "Regression end time          = ${REG_END_TIME}"
echo "Tip displacement difference  < ${TIP_TOL}"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    # Copy only the case inputs: results, logs and generated meshes from any
    # previous run of the tutorial are deliberately left behind
    for item in Allrun Allclean 0 constant system geometry \
        makeFluidSurface.py
    do
        cp -a "${SCRIPT_DIR}/${item}" "${CASE_DIR}/"
    done

    rm -rf "${CASE_DIR}/constant/polyMesh" \
           "${CASE_DIR}/constant/fluid/polyMesh" \
           "${CASE_DIR}/constant/solid/polyMesh" \
           "${CASE_DIR}/constant/triSurface"

    # Portable in-place edit (GNU and BSD sed)
    sed -i.orig \
        -e "s/^\(endTime[[:space:]]*\).*/\1${REG_END_TIME};/" \
        -e "s/^\(writeInterval[[:space:]]*\).*/\1${REG_END_TIME};/" \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.orig"
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

if [ "$CHECK_ONLY" = false ]; then
    prepare_case
    ( cd "${CASE_DIR}" && ./Allrun coarse > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allrun"
fi

# The case requires PETSc, cartesianMesh and python3: when any is missing
# Allrun writes a declared skip message, which is the only valid reason to
# skip the checks
if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

if [[ ! -f "${CASE_DIR}/${TIP_FILE}" ]]; then
    echo "FAIL: the case did not run in this environment:"
    echo "      expected output is missing and the tutorial did not declare a skip"
    echo "      (see ${CASE_DIR}/${ALLRUN_LOGFILE})"
    exit 1
fi

if ! grep -q "^End" "${CASE_DIR}/log.solids4Foam"; then
    echo "FAIL: solids4Foam did not finish (see ${CASE_DIR}/log.solids4Foam)"
    exit 1
fi

# Last line of the tip history: time Dx Dy Dz |D|
last=$(awk '($1 + 0) == $1 { line = $0 } END { print line }' \
    "${CASE_DIR}/${TIP_FILE}")
tip_time=$(echo "${last}" | awk '{print $1}')
tip_dy=$(echo "${last}" | awk '{print $3}')
tip_dz=$(echo "${last}" | awk '{print $4}')

if [[ -z "${tip_time}" ]] || \
    ! awk "BEGIN {exit !(${tip_time} + 0 >= ${REG_END_TIME} - 1e-9)}"
then
    echo "FAIL: the tip history stops at t = ${tip_time:-none}, short of the"
    echo "      requested end time ${REG_END_TIME}"
    exit 1
fi

failures=0
for pair in "Dy ${tip_dy} ${REF_TIP_DY}" "Dz ${tip_dz} ${REF_TIP_DZ}"; do
    IFS=' ' read -r name value ref <<< "${pair}"
    if awk -v v="${value}" -v r="${ref}" -v t="${TIP_TOL}" \
        'BEGIN { d = v - r; if (d < 0) d = -d; exit !(d < t && v == v) }'
    then
        printf "PASS: tip %s = %.6g (reference %.6g)\n" "${name}" "${value}" "${ref}"
    else
        printf "FAIL: tip %s = %.6g (reference %.6g)\n" "${name}" "${value}" "${ref}"
        failures=$((failures + 1))
    fi
done

if [ "${failures}" -ne 0 ]; then
    echo
    echo "Regression test FAILED (${failures} failure(s))"
    exit 1
fi

echo
echo "Regression test PASSED"
