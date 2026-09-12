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
# compressedSpheres regression test
# Runs the default frictionless case on the shipped coarse mesh up to a
# plate displacement of 2.5 mm, i.e. before the outer sphere buckles, and
# checks the compression force on the outer sphere (patch R_top):
#   1. against the value this test was calibrated with (regression check);
#   2. against the band spanned by the published solutions in
#      referenceData/frictionless, interpolated to 2.5 mm and widened by
#      REFERENCE_MARGIN (benchmark check).
# ============================================================

REGRESSION_END_TIME=0.25
REGRESSION_DISP_MM=2.5
EXPECTED_TIME_STEPS=25

# Calibrated with foam-extend-4.1: force_z = 0.999353 N at 2.5 mm
FORCE_Z_MIN=0.99
FORCE_Z_MAX=1.01

# Relative margin added to the reference band
REFERENCE_MARGIN=0.05

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
FORCE_FILE="postProcessing/0/solidForcesR_top.dat"
REFERENCE_DIR="referenceData/frictionless"

echo "============================================================"
echo "compressedSpheres regression test"
echo "Force_z on R_top at ${REGRESSION_DISP_MM} mm in [${FORCE_Z_MIN}, ${FORCE_Z_MAX}] N"
echo "and within the published solutions +/- $(awk "BEGIN {print 100*${REFERENCE_MARGIN}}")%"
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

    sed -i.bak "s/^endTime[[:space:]]\+0.75;/endTime         ${REGRESSION_END_TIME};/" \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"
}

# Linearly interpolate a digitised reference curve at displacement x0 (mm)
# Arguments: file, x0, x scale factor to mm, y scale factor to N
interpolate_reference() {
    awk -v x0="$2" -v sx="$3" -v sy="$4" '
        !/^#/ && NF >= 2 {
            x = $1*sx; y = $2*sy
            if (n && px <= x0 && x >= x0 && x > px) {
                print py + (y - py)*(x0 - px)/(x - px); found = 1; exit
            }
            px = x; py = y; n++
        }
        END { if (!found) exit 1 }' "$1"
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

failures=0

# Solver log checks
if [[ ! -f "${CASE_DIR}/${SOLVER_LOGFILE}" ]]; then
    echo "FAIL: Could not find ${SOLVER_LOGFILE}"
    exit 1
fi

if grep -Eq "FOAM FATAL|^ERROR$|\[stack trace\]" "${CASE_DIR}/${SOLVER_LOGFILE}"; then
    echo "FAIL: ${SOLVER_LOGFILE} reports an error"
    failures=$((failures + 1))
elif ! grep -q "^End" "${CASE_DIR}/${SOLVER_LOGFILE}"; then
    echo "FAIL: ${SOLVER_LOGFILE} does not finish with End"
    failures=$((failures + 1))
else
    echo "PASS: solids4Foam finished without errors"
fi

if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
    echo "FAIL: Could not find ${FORCE_FILE}"
    exit 1
fi

# Only the entries after time zero are counted
n_steps=$(awk '!/^#/ && NF && $1 + 0 != 0 {count++} END {print count + 0}' \
    "${CASE_DIR}/${FORCE_FILE}")
final_time=$(awk '!/^#/ && NF {time = $1} END {print time}' \
    "${CASE_DIR}/${FORCE_FILE}")
final_force_z=$(awk '!/^#/ && NF {force = $4} END {print force}' \
    "${CASE_DIR}/${FORCE_FILE}")

if [[ "${n_steps}" != "${EXPECTED_TIME_STEPS}" ]] \
|| ! awk "BEGIN {exit !((${final_time} - ${REGRESSION_END_TIME})^2 < 1e-12)}"; then
    echo "FAIL: ${n_steps} force entries up to time ${final_time}; expected ${EXPECTED_TIME_STEPS} up to ${REGRESSION_END_TIME}"
    failures=$((failures + 1))
fi

if [[ -z "${final_force_z}" ]]; then
    echo "FAIL: Could not extract the final force_z"
    exit 1
fi

# 1. Regression check
if awk "BEGIN {exit !(${final_force_z} >= ${FORCE_Z_MIN} && ${final_force_z} <= ${FORCE_Z_MAX})}"; then
    printf "PASS: Final force_z = %.6g N\n" "${final_force_z}"
else
    printf "FAIL: Final force_z = %.6g N; expected between %g and %g N\n" \
        "${final_force_z}" "${FORCE_Z_MIN}" "${FORCE_Z_MAX}"
    failures=$((failures + 1))
fi

# 2. Benchmark check against the published solutions
# The unit conversions are given in the header of each reference data file
REF_DIR="${SCRIPT_DIR}/${REFERENCE_DIR}"
ref_values=(
    "$(interpolate_reference "${REF_DIR}/Abaqus.dat" "${REGRESSION_DISP_MM}" 1 1)"
    "$(interpolate_reference "${REF_DIR}/FEBio.dat" "${REGRESSION_DISP_MM}" 10 -1)"
    "$(interpolate_reference "${REF_DIR}/Areias.dat" "${REGRESSION_DISP_MM}" 10 -1)"
    "$(interpolate_reference "${REF_DIR}/PusoLaursen.dat" "${REGRESSION_DISP_MM}" 10 -1)"
)

ref_min=$(printf "%s\n" "${ref_values[@]}" | sort -g | head -1)
ref_max=$(printf "%s\n" "${ref_values[@]}" | sort -g | tail -1)
band_min=$(awk "BEGIN {print ${ref_min}*(1 - ${REFERENCE_MARGIN})}")
band_max=$(awk "BEGIN {print ${ref_max}*(1 + ${REFERENCE_MARGIN})}")

printf "Published solutions at %s mm (Abaqus, FEBio, Areias, Puso-Laursen): %s N\n" \
    "${REGRESSION_DISP_MM}" "$(printf "%.4g " "${ref_values[@]}" | sed 's/ $//')"

if awk "BEGIN {exit !(${final_force_z} >= ${band_min} && ${final_force_z} <= ${band_max})}"; then
    printf "PASS: Final force_z = %.6g N lies within the reference band [%.4g, %.4g] N\n" \
        "${final_force_z}" "${band_min}" "${band_max}"
else
    printf "FAIL: Final force_z = %.6g N lies outside the reference band [%.4g, %.4g] N\n" \
        "${final_force_z}" "${band_min}" "${band_max}"
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
