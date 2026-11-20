#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"
source "${SOLIDS4FOAM_SCRIPTS}"

# GNU sed, for the in-place edit below
solids4Foam::requireGnuSed

# The end-to-end check of the non-local Lemaitre damage law: the damage
# regularised by the implicit-gradient field damageNonLocal, which the
# updated Lagrangian solid model solves, against the Abaqus reaction force of
# Whelan et al. (2025)
# ============================================================
# notchedBarLemaitre regression test
# Checks the loading force at a displacement of 0.6 mm, past the peak, where
# the damage softens the bar.
# ============================================================

# The force at 0.6 mm. Abaqus gives 27.18 kN, and this mesh 1.8% less; the
# local form of the law gives 16.4 kN
END_TIME=0.6
ABAQUS_FORCE=27.18
ABAQUS_REL_TOL=3e-2

# OpenFOAM.com v2512's, at END_TIME
REF_FORCE=26.676648
REF_REL_TOL=1e-3

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "notchedBarLemaitre regression test"
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

if grep -q "nonLocal = classic" "${CASE_DIR}/log.solids4Foam" 2>/dev/null \
 && grep -q "Solving for damageNonLocal" "${CASE_DIR}/log.solids4Foam"
then
    echo "PASS: the law is non-local and damageNonLocal is solved"
else
    echo "FAIL: the non-local law or the damageNonLocal equation is missing"
    failures=$((failures + 1))
fi

# The non-local damage at END_TIME: written, and damaging without failing
damage_file="${CASE_DIR}/${END_TIME}/damageNonLocal"
max_damage=""
if [[ -f "${damage_file}" ]]; then
    max_damage=$(awk '
        /^internalField/ {inField = 1; next}
        inField && /^\(/ {inList = 1; next}
        inList && /^\)/ {exit}
        inList {if ($1 > m) m = $1}
        END {printf "%.6g", m}
    ' "${damage_file}")
fi

if [[ -n "${max_damage}" ]] \
 && awk "BEGIN {exit !(${max_damage} > 0.1 && ${max_damage} < 0.99)}"
then
    printf "PASS: maximum non-local damage %.4g\n" "${max_damage}"
else
    printf "FAIL: maximum non-local damage '%s' not in (0.1, 0.99)\n" \
        "${max_damage}"
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

    if awk "BEGIN {exit !(${final_disp} > 0.599 && ${final_disp} < 0.601)}"
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
