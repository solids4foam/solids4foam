#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"
source "${SOLIDS4FOAM_SCRIPTS}"

# GNU sed, for the in-place edit below
solids4Foam::requireGnuSed

# The end-to-end check of the logStrain law: a small-strain J2 law driven by
# the elastic logarithmic strain, under the updated Lagrangian solid model,
# against the Abaqus reaction force of Whelan et al. (2025)
# ============================================================
# notchedBar regression test
# Checks the loading force at a displacement of 0.3 mm.
# ============================================================

# The force at 0.3 mm. Abaqus gives 32.746 kN; the logarithmic strain is
# expected slightly below it, as Whelan et al. report, and this coarse mesh
# gives 0.4% less
END_TIME=0.3
ABAQUS_FORCE=32.746
ABAQUS_REL_TOL=1.5e-2

# OpenFOAM.com v2512's, at END_TIME
REF_FORCE=32.623992
REF_REL_TOL=5e-4

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "notchedBar regression test"
echo "Loading force at ${END_TIME} mm within ${ABAQUS_REL_TOL} of Abaqus"
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

    "${SOLIDS4FOAM_SED}" -i.bak "s/^endTime         1;/endTime         ${END_TIME};/" \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"

    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if grep -q "Logarithmic strain on top of linearElasticMisesPlastic" \
    "${CASE_DIR}/log.solids4Foam" 2>/dev/null
then
    echo "PASS: took its material from the logStrain law"
else
    echo "FAIL: did not construct the logStrain law"
    failures=$((failures + 1))
fi

force_file=$(find "${CASE_DIR}/postProcessing" \
    -name 'solidForcesDisplacementsup.dat' -print 2> /dev/null | tail -n 1)

final_force=""
if [[ -z "${force_file}" ]]; then
    echo "FAIL: could not find the loading force history"
    failures=$((failures + 1))
else
    # The force on the 5 degree wedge, scaled to the whole bar, in kN
    final_force=$(awk 'END {printf "%.10g", $6*72e-3}' "${force_file}")
    final_disp=$(awk 'END {print $3*1e3}' "${force_file}")

    if awk "BEGIN {exit !(${final_disp} > 0.299 && ${final_disp} < 0.301)}"
    then
        printf "PASS: reached a displacement of %.4g mm\n" "${final_disp}"
    else
        printf "FAIL: stopped at a displacement of %.4g mm\n" "${final_disp}"
        failures=$((failures + 1))
    fi

    a="${ABAQUS_FORCE}"; b="${final_force}"
    if awk "BEGIN {exit !(($a - $b)^2 <= (${ABAQUS_REL_TOL}*$a)^2)}"; then
        printf "PASS: force agrees with Abaqus (%.6g vs %.6g kN)\n" "$b" "$a"
    else
        printf "FAIL: force differs from Abaqus (%.6g vs %.6g kN)\n" "$b" "$a"
        failures=$((failures + 1))
    fi

    a="${REF_FORCE}"
    if awk "BEGIN {exit !(($a - $b)^2 <= (${REF_REL_TOL}*$a)^2)}"; then
        printf "PASS: force matches the reference (%.8g vs %.8g kN)\n" "$b" "$a"
    else
        printf "FAIL: force differs from the reference (%.8g vs %.8g kN)\n" "$b" "$a"
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
