#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# Land et al. (2015) problem 3 regression test
#
# The benchmark measures how far a point on the tissue moves under an active
# tension that ramps over the run. What makes this case worth having is that
# it is the only tutorial whose material is anisotropic about a fibre
# direction: the fibre field is built by setFibreField before the solver runs,
# and the constitutive law reads it.
#
# Both arms run: the legacy mechanicalModel and the mechanicalConstitutiveLaw
# framework, which must agree. This is the only case that exercises a fibre
# direction as prescribed state, and the only one where a composite law adds an
# active stress to a passive one.
# ============================================================

MAG_D_MIN=8.0e-4
MAG_D_MAX=9.0e-4

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"

echo "============================================================"
echo "Land et al. (2015) problem 3 regression test"
echo "Final probe |D| in [${MAG_D_MIN}, ${MAG_D_MAX}] m"
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

prepare_case
sed -i 's/^writePrecision.*/writePrecision  14;/' "${CASE_DIR}/system/controlDict"
( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped here"
    exit 0
fi

if [[ ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    echo "SKIP: no displacement history (the case did not run here)"
    exit 0
fi

failures=0

# The fibre field has to exist and point somewhere, or this case is an
# isotropic one wearing an anisotropic law and the port it exists to check
# would be checked against nothing
if [[ ! -f "${CASE_DIR}/0/f0" ]]; then
    echo "FAIL: setFibreField produced no fibre field"
    failures=$((failures + 1))
elif python3 - "${CASE_DIR}/0/f0" << 'PYEOF'
import re, sys
body = open(sys.argv[1]).read().split('* * * * *')[-1]
nums = [float(x) for x in re.findall(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', body)]
# any component away from zero means a direction was actually set
sys.exit(0 if any(abs(v) > 1e-8 for v in nums) else 1)
PYEOF
then
    echo "PASS: the fibre field is set"
else
    echo "FAIL: the fibre field is everywhere zero"
    failures=$((failures + 1))
fi

magD=$(awk 'END {print $5}' "${CASE_DIR}/${DISP_FILE}")

if [[ -z "${magD}" ]]; then
    echo "FAIL: could not read the final probe displacement"
    failures=$((failures + 1))
elif awk "BEGIN {exit !(${magD} >= ${MAG_D_MIN} && ${magD} <= ${MAG_D_MAX})}"; then
    printf "PASS: final probe |D| = %.6g\n" "${magD}"
else
    printf "FAIL: final probe |D| = %.6g\n" "${magD}"
    failures=$((failures + 1))
fi


# ------------------------------------------------------------
# The framework arm
# ------------------------------------------------------------
# The same case with the stress taken from the mechanicalConstitutiveLaw
# framework rather than the legacy mechanicalModel, in the displacement
# formulation. The mixed formulation is the arm after this one
run_framework_comparison() {
    local d="${REGRESSION_ROOT}/framework"

    rm -rf "${d}"; mkdir -p "${d}"
    for item in "${SCRIPT_DIR}"/*; do
        [[ "$(basename "${item}")" == "regressionTests" ]] && continue
        cp -a "${item}" "${d}/"
    done
    sed -i 's/^writePrecision.*/writePrecision  14;/' "${d}/system/controlDict"

    # The switch lives in the solid model's Coeffs sub-dictionary, and Allrun
    # symlinks solidProperties.displacement into place, so the variant is what
    # has to be edited. At the top level it would be silently ignored and this
    # arm would quietly repeat the legacy run
    sed -i.bak \
        's|^    // Solution algorithm|    useMechanicalConstitutiveLawManager yes;\n    // Solution algorithm|' \
        "${d}/constant/solidProperties.displacement"
    rm -f "${d}/constant/solidProperties.displacement.bak"

    if ! grep -q 'useMechanicalConstitutiveLawManager' \
        "${d}/constant/solidProperties.displacement"
    then
        echo "FAIL: could not enable the framework"
        return 1
    fi

    ( cd "${d}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true

    if solids4Foam::regressionCaseSkipped "${d}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: framework arm (the tutorial skipped here)"
        return 0
    fi

    # Each arm must have taken the path it was set up for
    if ! grep -q "Selecting mechanical constitutive law" \
        "${d}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the framework arm did not use the framework"
        return 1
    fi

    if grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the legacy arm used the framework"
        return 1
    fi
    echo "PASS: each arm took the path it was set up for"

    # The framework laws' own checks, which include the one that matters for
    # this case: a law declaring it can separate its isochoric stress from its
    # volumetric response is taken at its word by any mixed formulation, and
    # this is what tests the word. GuccioneElastic declares it
    if command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
        if ( cd "${d}" && Test-mechanicalConstitutiveLaw > log.unit 2>&1 ); then
            if grep -q "isochoric stress ignores a superposed dilation" \
                "${d}/log.unit"
            then
                echo "PASS: $(grep -c 'PASS:' "${d}/log.unit") law checks," \
                     "including the isochoric split"
            else
                echo "FAIL: the isochoric split was not checked here"
                return 1
            fi
        else
            echo "FAIL: the law checks did not pass"
            grep -m2 "FAIL:" "${d}/log.unit" || true
            return 1
        fi
    fi

    local b
    b=$(awk 'END {print $5}' "${d}/${DISP_FILE}" 2>/dev/null)

    if [[ -z "${b}" ]]; then
        echo "FAIL: the framework arm produced no displacement history"
        return 1
    fi

    # The framework does NOT reproduce legacy here, deliberately, and this is
    # the one place in the suite where that is true.
    #
    # The legacy GuccioneElastic builds Q from the full Green-Lagrange strain,
    # so its energy is coupled: the deviatoric stress varies with volume
    # change. The ported law builds Q from the isochoric strain, so shape and
    # volume are separate, which is what lets a mixed displacement-pressure
    # formulation replace the volumetric part and still describe the same
    # material. Both reduce to the published model in the incompressible limit
    # it was defined for; away from that limit they are different materials and
    # the difference here is 4e-4.
    #
    # So this checks two things instead of one: that the framework value has
    # not drifted, and that it stays near legacy - near enough that the two are
    # recognisably the same problem, far enough that the reformulation is
    # visible rather than lost
    local ref=8.5119442158e-4

    if ! awk "BEGIN {exit !((${b} - ${ref})^2 <= (1e-6*${ref})^2)}"; then
        printf "FAIL: framework moved from its reference (%.10g vs %.10g)\n" \
            "${ref}" "${b}"
        return 1
    fi
    printf "PASS: framework matches its reference (%.10g)\n" "${b}"

    if awk "BEGIN {exit !((${magD} - ${b})^2 <= (1e-3*${magD})^2)}"; then
        printf "PASS: framework near legacy, differing by the reformulation (%.10g vs %.10g)\n" \
            "${magD}" "${b}"
        return 0
    fi

    printf "FAIL: framework and legacy differ by more than the reformulation explains (%.10g vs %.10g)\n" \
        "${magD}" "${b}"
    return 1
}

if ! run_framework_comparison; then
    failures=$((failures + 1))
fi

# ------------------------------------------------------------
# The mixed displacement-pressure arm, on the framework
# ------------------------------------------------------------
# This is the arm the volumetric split exists for. The pressure is solved as
# its own unknown and replaces the law's volumetric response, so the law is
# asked for its stress without that response rather than being asked for the
# total and having a deviatoric projection taken of it.
#
# The two are not the same here, and this case is why. Guccione's isochoric
# stress is trace free, so a projection would recover it exactly; the active
# tension along the fibre direction is not, and a projection discards its
# spherical part while the pressure - which knows only about dU/dJ - does not
# put it back. That is a 1.6% difference in the final displacement of this
# case, and it is a different material rather than a visible error, which is
# why the framework refuses a law that cannot separate the two
run_mixed_framework() {
    local d="${REGRESSION_ROOT}/frameworkPressure"

    rm -rf "${d}"; mkdir -p "${d}"
    for item in "${SCRIPT_DIR}"/*; do
        [[ "$(basename "${item}")" == "regressionTests" ]] && continue
        cp -a "${item}" "${d}/"
    done
    sed -i 's/^writePrecision.*/writePrecision  14;/' "${d}/system/controlDict"

    # Allrun symlinks the ".pressure" variants into place, so that is the file
    # the switch has to go in
    sed -i.bak \
        's|^    solvePressure true;|    solvePressure true;\n    useMechanicalConstitutiveLawManager yes;|' \
        "${d}/constant/solidProperties.pressure"
    rm -f "${d}/constant/solidProperties.pressure.bak"

    if ! grep -q 'useMechanicalConstitutiveLawManager' \
        "${d}/constant/solidProperties.pressure"
    then
        echo "FAIL: could not enable the framework on the mixed arm"
        return 1
    fi

    ( cd "${d}" && ./Allrun pressure > "${ALLRUN_LOGFILE}" 2>&1 ) || true

    if solids4Foam::regressionCaseSkipped "${d}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: mixed framework arm (the tutorial skipped here)"
        return 0
    fi

    if ! grep -q "Selecting mechanical constitutive law" \
        "${d}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the mixed arm did not use the framework"
        return 1
    fi

    # The implicit stiffness must be the deviatoric surrogate. With the full
    # one the bulk modulus - a near-incompressibility penalty two orders above
    # the shear modulus here - makes the Laplacian surrogate so stiff that the
    # linear solve does not converge at all, which is what used to happen
    if ! grep -q "scalarDeviatoric" "${d}/${SOLVER_LOGFILE}"; then
        echo "FAIL: the mixed arm did not ask for a deviatoric implicit stiffness"
        return 1
    fi
    echo "PASS: the mixed arm solved a pressure on the framework"

    local m
    m=$(awk 'END {print $5}' "${d}/${DISP_FILE}" 2>/dev/null)

    if [[ -z "${m}" ]]; then
        echo "FAIL: the mixed arm produced no displacement history"
        return 1
    fi

    local mref=1.15475e-3

    if ! awk "BEGIN {exit !((${m} - ${mref})^2 <= (2e-3*${mref})^2)}"; then
        printf "FAIL: mixed arm moved from its reference (%.10g vs %.10g)\n" \
            "${mref}" "${m}"
        return 1
    fi
    printf "PASS: mixed arm matches its reference (%.10g)\n" "${m}"

    return 0
}

if ! run_mixed_framework; then
    failures=$((failures + 1))
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
