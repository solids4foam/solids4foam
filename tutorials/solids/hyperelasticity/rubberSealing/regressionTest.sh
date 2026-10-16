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
# rubberSealing regression test
# Checks the vertical force on the compressed (top) surface and the
# displacement magnitude of the outer upper corner of the inclined wall at the
# end of the loading, plus a loose sanity check of the peak equivalent stress.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

# Self-regression values, not validated reference values: solids4foam gives a
# final force_y of -0.171263 N and a final probe displacement magnitude of
# 0.000975391 m (foam-extend-4.1); the bands are +/- 1% about these values
FORCE_Y_MIN=-0.17298
FORCE_Y_MAX=-0.16955
DISP_MAG_MIN=0.00096564
DISP_MAG_MAX=0.00098515

# Loose sanity check only: the final peak equivalent stress (903.8 kPa with
# foam-extend-4.1) must lie within 10% of the 882 kPa reported by Pascon
# (2019). This is not a validation, as Pascon uses an incompressible Yeoh-type
# model with a different initial shear modulus and a finite-strain plane
# stress formulation, whereas this case uses a neo-Hookean law with the
# solids4foam (linear) plane stress approximation
SIGMA_EQ_MIN=793800
SIGMA_EQ_MAX=970200

# The displacement is fully applied at the end time
END_TIME=100

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
FORCE_FILE="postProcessing/0/solidForcestop.dat"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"

echo "============================================================"
echo "rubberSealing regression test"
echo "Final top force_y in [${FORCE_Y_MIN}, ${FORCE_Y_MAX}] N"
echo "Final probe disp magnitude in [${DISP_MAG_MIN}, ${DISP_MAG_MAX}] m"
echo "Final peak sigmaEq in [${SIGMA_EQ_MIN}, ${SIGMA_EQ_MAX}] Pa (loose check)"
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

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

if [[ ! -f "${CASE_DIR}/${SOLVER_LOGFILE}" ]]; then
    echo "FAIL: Could not find ${SOLVER_LOGFILE}"
    exit 1
fi

if grep -qE 'FOAM FATAL|^ERROR$|\[stack trace\]' "${CASE_DIR}/${SOLVER_LOGFILE}" \
    || ! grep -q '^End' "${CASE_DIR}/${SOLVER_LOGFILE}"
then
    echo "FAIL: solids4Foam did not complete"
    exit 1
fi

# The solver may print End even if the solution diverged, so also check the
# log, the written fields and the histories for nan
if grep -qiw 'nan' "${CASE_DIR}/${SOLVER_LOGFILE}" \
    || grep -rqiw 'nan' "${CASE_DIR}"/[0-9]* "${CASE_DIR}/postProcessing"
then
    echo "FAIL: nan found in the solver log, fields or histories"
    exit 1
fi

if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" || ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    echo "FAIL: Could not find one or more history files"
    exit 1
fi

final_time=$(awk 'END {print $1}' "${CASE_DIR}/${FORCE_FILE}")
final_force_y=$(awk 'END {print $3}' "${CASE_DIR}/${FORCE_FILE}")
final_disp_mag=$(awk 'END {print $5}' "${CASE_DIR}/${DISP_FILE}")
final_sigma_eq=$(awk '/^Max sigmaEq/ {v = $NF} END {print v}' "${CASE_DIR}/${SOLVER_LOGFILE}")

if [[ -z "${final_force_y}" || -z "${final_disp_mag}" || -z "${final_sigma_eq}" ]]; then
    echo "FAIL: Could not extract final force/displacement/equivalent stress"
    exit 1
fi

failures=0

if ! awk "BEGIN {exit !(${final_time} == ${END_TIME})}"; then
    echo "FAIL: Final time = ${final_time}, expected ${END_TIME}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_force_y} >= ${FORCE_Y_MIN} && ${final_force_y} <= ${FORCE_Y_MAX})}"; then
    printf "PASS: Final force_y = %.6g\n" "${final_force_y}"
else
    printf "FAIL: Final force_y = %.6g\n" "${final_force_y}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_disp_mag} >= ${DISP_MAG_MIN} && ${final_disp_mag} <= ${DISP_MAG_MAX})}"; then
    printf "PASS: Final disp magnitude = %.6g\n" "${final_disp_mag}"
else
    printf "FAIL: Final disp magnitude = %.6g\n" "${final_disp_mag}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_sigma_eq} >= ${SIGMA_EQ_MIN} && ${final_sigma_eq} <= ${SIGMA_EQ_MAX})}"; then
    printf "PASS: Final peak sigmaEq = %.6g (loose check)\n" "${final_sigma_eq}"
else
    printf "FAIL: Final peak sigmaEq = %.6g (loose check)\n" "${final_sigma_eq}"
    failures=$((failures + 1))
fi

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
