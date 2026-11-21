#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"
source "${SOLIDS4FOAM_SCRIPTS}"

# GNU sed, for the in-place edit below
solids4Foam::requireGnuSed

# The end-to-end check of the phase field fracture law under logStrain, with
# the phase field solved by the updated Lagrangian solid model as an
# implicit-gradient field with a coefficient, against the normalised stress
# of Borden et al. (2016)
# ============================================================
# flatNotchedBarPhaseField regression test
# Checks the normalised stress at a normalised displacement of 0.0042, at
# the peak, where the phase field has grown but the crack has not run.
# ============================================================

# The normalised stress at 0.0042. Borden et al. give 348.2 MPa, and this
# mesh 0.1% more
END_TIME=0.28
BORDEN_STRESS=348.2
BORDEN_REL_TOL=2e-2

# OpenFOAM.com v2512's, at END_TIME
REF_STRESS=348.66581
REF_REL_TOL=1e-3

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "flatNotchedBarPhaseField regression test"
echo "Normalised stress at 0.0042 within ${BORDEN_REL_TOL} of Borden et al."
echo "============================================================"
echo

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

failures=0
CASE_DIR="${REGRESSION_ROOT}/main"

if [ "$CHECK_ONLY" = false ]; then
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"
    for item in "${SCRIPT_DIR}"/*; do
        if [[ "$(basename "${item}")" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done

    "${SOLIDS4FOAM_SED}" -i.bak "s/^endTime         0.32;/endTime         ${END_TIME};/" \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"

    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if grep -q "Phase field fracture: Gc = 60000" "${CASE_DIR}/log.solids4Foam" \
    2>/dev/null \
 && grep -q "Solving for phaseField" "${CASE_DIR}/log.solids4Foam"
then
    echo "PASS: the phase field law is read and phaseField is solved"
else
    echo "FAIL: the phase field law or the phaseField equation is missing"
    failures=$((failures + 1))
fi

# The phase field at END_TIME: written, and grown without a crack
phase_file="${CASE_DIR}/${END_TIME}/phaseField"
max_phase=""
if [[ -f "${phase_file}" ]]; then
    max_phase=$(awk '
        /^internalField/ {inField = 1; next}
        inField && /^\(/ {inList = 1; next}
        inList && /^\)/ {exit}
        inList {if ($1 > m) m = $1}
        END {printf "%.6g", m}
    ' "${phase_file}")
fi

if [[ -n "${max_phase}" ]] \
 && awk "BEGIN {exit !(${max_phase} > 0.01 && ${max_phase} < 0.5)}"
then
    printf "PASS: maximum phase field %.4g\n" "${max_phase}"
else
    printf "FAIL: maximum phase field '%s' not in (0.01, 0.5)\n" \
        "${max_phase}"
    failures=$((failures + 1))
fi

force_file=$(find "${CASE_DIR}/postProcessing" \
    -name 'solidForcesDisplacementsup.dat' -print 2> /dev/null | tail -n 1)

final_stress=""
if [[ -z "${force_file}" ]]; then
    echo "FAIL: could not find the loading force history"
    failures=$((failures + 1))
else
    # The force on the quarter plate over the area of its notched section,
    # 20.4768 mm^2, in MPa, against the displacement over 76.2 mm
    final_stress=$(awk 'END {printf "%.10g", $6/20.4768}' "${force_file}")
    final_disp=$(awk 'END {printf "%.6g", $3/76.2e-3}' "${force_file}")

    if awk "BEGIN {exit !(${final_disp} > 0.004199 && ${final_disp} < 0.004201)}"
    then
        printf "PASS: reached a normalised displacement of %.6g\n" \
            "${final_disp}"
    else
        printf "FAIL: stopped at a normalised displacement of %.6g\n" \
            "${final_disp}"
        failures=$((failures + 1))
    fi

    a="${BORDEN_STRESS}"; b="${final_stress}"
    if awk "BEGIN {exit !(($a - $b)^2 <= (${BORDEN_REL_TOL}*$a)^2)}"; then
        printf "PASS: stress agrees with Borden et al. (%.6g vs %.6g MPa)\n" \
            "$b" "$a"
    else
        printf "FAIL: stress differs from Borden et al. (%.6g vs %.6g MPa)\n" \
            "$b" "$a"
        failures=$((failures + 1))
    fi

    a="${REF_STRESS}"
    if awk "BEGIN {exit !(($a - $b)^2 <= (${REF_REL_TOL}*$a)^2)}"; then
        printf "PASS: stress matches the reference (%.8g vs %.8g MPa)\n" \
            "$b" "$a"
    else
        printf "FAIL: stress differs from the reference (%.8g vs %.8g MPa)\n" \
            "$b" "$a"
        failures=$((failures + 1))
    fi
fi

if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
fi

echo "============================================================"
if [[ ${failures} -eq 0 ]]; then
    echo "Regression test PASSED"
    echo "============================================================"
    exit 0
else
    echo "Regression test FAILED (${failures} checks)"
    echo "============================================================"
    exit 1
fi
