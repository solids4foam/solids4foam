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
# beamCantilever regression test
# Checks the free-end displacement against the analytical solution for
# a cantilever bent by a prescribed end rotation, and against beamFoam.
# ============================================================

# Beam length (m) and final prescribed end rotation (rad)
L=10
PSI_END=6.283185307179586

# Largest allowed distance from the analytical free-end position, relative to
# L, over the whole history; the 50-segment mesh gives about 2.4e-4
ANALYTICAL_TOL=5e-4

# The reference values are the standalone beamFoam answer for this case
# (beamFoam 9263c5c, the commit in src/beamModels/BEAMFOAM_VERSION), which
# solids4Foam reproduced in every written digit on OpenFOAM.com v2412. The
# history is written to six figures, so the tolerance, 1e-5 relative, is a few
# units in the last digit
REF_TIME=0.75
REF_WX=-12.1205
REF_WZ=2.12051
REF_REL_TOL=1e-5

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "beamCantilever regression test"
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
# Against the analytical solution
#
# An end rotation psi bends the beam into a circular arc. The free end moves
# by w_x = -(L - L sin(psi)/psi) along the beam and w_z = L (1 - cos(psi))/psi
# across it; the rotation increases linearly to PSI_END at t = 1
# ------------------------------------------------------------

max_err=$(awk -v L="${L}" -v psiEnd="${PSI_END}" '
    NR > 1 && $1 > 0 {
        p = psiEnd*$1
        ax = -(L - L*sin(p)/p)
        az = L*(1 - cos(p))/p
        e = sqrt(($2 - ax)^2 + ($4 - az)^2)/L
        if (e > max) max = e
    }
    END {print max + 0}
' "${history_file}")

if awk "BEGIN {exit !(${max_err} <= ${ANALYTICAL_TOL})}"; then
    printf "PASS: max distance from the analytical arc = %.3g L\n" "${max_err}"
else
    printf "FAIL: max distance from the analytical arc = %.3g L (allowed %g L)\n" \
        "${max_err}" "${ANALYTICAL_TOL}"
    failures=$((failures + 1))
fi

# ------------------------------------------------------------
# Against beamFoam
# ------------------------------------------------------------

IFS=" " read -r wx wy wz wmag <<< "$(tip_at "${REF_TIME}")"
if [[ -z "${wx:-}" ]]; then
    echo "FAIL: no displacement written at t = ${REF_TIME}"
    exit 1
fi

check_value "w_x(${REF_TIME})" "${wx}" "${REF_WX}" "${L}" "${REF_REL_TOL}"
check_value "w_z(${REF_TIME})" "${wz}" "${REF_WZ}" "${L}" "${REF_REL_TOL}"

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
