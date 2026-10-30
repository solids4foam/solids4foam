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
# beamFollowerForce regression test
# Checks the free-end displacement against the Simo and Vu-Quoc (1986)
# reference solution and against beamFoam.
# ============================================================

# Beam length (m) and final follower load (kN)
L=100
LOAD_END=134

# The last point of the Simo and Vu-Quoc solution in referenceResultsSimo.dat:
# w_z (m) at 129.73 kN. beamFoam is within about 1.4% of L there; allow 3%
SIMO_LOAD=129.7349161
SIMO_WZ=-48.82606369
SIMO_TOL=0.03

# The reference values are the standalone beamFoam answer for this case
# (beamFoam 9263c5c, the commit in src/beamModels/BEAMFOAM_VERSION), which
# solids4Foam reproduced in every written digit on OpenFOAM.com v2412. The
# history is written to six figures, so the tolerance, 1e-5 relative, is a few
# units in the last digit
REF_TIME=1
REF_WX=-79.7138
REF_WZ=-49.3637
REF_MAG=93.7606
REF_REL_TOL=1e-5

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "beamFollowerForce regression test"
echo "============================================================"
echo

prepare_case() {
    local d="$1"
    rm -rf "${d}"
    mkdir -p "${d}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${d}/"
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
    prepare_case "${CASE_DIR}"
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

history_file="${CASE_DIR}/postProcessing/0/beamDisplacements_right.dat"
if [[ ! -f "${history_file}" ]]; then
    echo "FAIL: Could not find ${history_file#"${CASE_DIR}"/}"
    exit 1
fi

failures=0

# Free-end displacement at a time: prints "w_x w_y w_z mag(w)"
tip_at() {
    awk -v t="$1" 'NR > 1 && $1 == t {print $2, $3, $4, $5}' "${history_file}"
}

# Compare one value with its reference, to a tolerance relative to a scale
check_value() {
    local name="$1" value="$2" ref="$3" scale="$4" tol="$5"

    if awk "BEGIN {d = ${value} - (${ref}); if (d < 0) d = -d;
                   exit !(d <= ${tol} * ${scale})}"
    then
        printf "PASS: %s = %.8g (reference %.8g)\n" "${name}" "${value}" "${ref}"
    else
        printf "FAIL: %s = %.8g (reference %.8g, tolerance %g x %g)\n" \
            "${name}" "${value}" "${ref}" "${tol}" "${scale}"
        failures=$((failures + 1))
    fi
}

# ------------------------------------------------------------
# Against Simo and Vu-Quoc (1986)
#
# The load increases linearly to LOAD_END at t = 1, so interpolate the
# history at the time of the reference load
# ------------------------------------------------------------

wz_simo=$(awk -v t="$(awk "BEGIN {print ${SIMO_LOAD}/${LOAD_END}}")" '
    NR > 1 {
        if ($1 >= t) {
            if (prevT == "") {print $4} else {print prevW + ($4 - prevW)*(t - prevT)/($1 - prevT)}
            exit
        }
        prevT = $1; prevW = $4
    }
' "${history_file}")

check_value "w_z at ${SIMO_LOAD} kN" "${wz_simo}" "${SIMO_WZ}" "${L}" "${SIMO_TOL}"

# ------------------------------------------------------------
# Against beamFoam
# ------------------------------------------------------------

IFS=" " read -r wx wy wz wmag <<< "$(tip_at "${REF_TIME}")"
if [[ -z "${wx:-}" ]]; then
    echo "FAIL: no displacement written at t = ${REF_TIME}"
    exit 1
fi

check_value "w_x(${REF_TIME})" "${wx}" "${REF_WX}" "${REF_MAG}" "${REF_REL_TOL}"
check_value "w_z(${REF_TIME})" "${wz}" "${REF_WZ}" "${REF_MAG}" "${REF_REL_TOL}"
check_value "|w|(${REF_TIME})" "${wmag}" "${REF_MAG}" "${REF_MAG}" "${REF_REL_TOL}"

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
