#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
DISP_FILE="postProcessing/0/solidPointDisplacement_pointHistory.dat"

SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
elif command -v solids4FoamScripts.sh > /dev/null 2>&1; then
    source solids4FoamScripts.sh
fi

if ! declare -F solids4Foam::regressionCaseSkipped > /dev/null 2>&1; then
    solids4Foam::regressionCaseSkipped() {
        local LOG_FILE="$1"
        [[ -f "${LOG_FILE}" ]] || return 1
        grep -Eq \
"This case currently only runs in foam-extend|\
This case currently does not run with foam-extend|\
This case currently does not run with OpenFOAM.org|\
Skipping this case as it does not currently working with OpenFOAM.org|\
OpenFOAM-v[0-9]+ or a newer version is required|\
Skipping this case as PETSc is not installed" \
            "${LOG_FILE}"
    }
fi

# ============================================================
# ratCarotid regression test
#
# An artery wall, two symmetric fibre families, inflated to 25 kPa.
#
# Two formulations, and only one of them works.
#
# The legacy one - coupledPressureDisplacementSolid with the legacy
# HolzapfelGasserOgdenElastic - is foam-extend only and does not reach the end
# time: it stalls at a relative residual of about 0.99 and dies at t = 0.68.
# That is longstanding rather than a regression; the same failure happens on
# the development branch this work started from. It is not run here.
#
# The framework one - nonLinearGeometryTotalLagrangianTotalDisplacement with
# solvePressure, taking its stress from the mechanicalConstitutiveLaw
# framework - runs to completion on any fork, which is what this checks.
#
# The two do not agree closely, and this test does not pretend they do. Over
# the range where the legacy run survives they differ by 2 to 4 per cent at
# high load and by up to 25 per cent at moderate load. That difference is not
# the material: the bulk modulus is a penalty here where the legacy law is
# exactly incompressible, and raising it a hundredfold moves the answer by
# 0.1 per cent, so this is at the incompressible limit already. Nor is it the
# pressure stabilisation, which moves it by 4 per cent over a sixteenfold
# sweep. What remains is the difference between an incremental updated
# Lagrangian solver on a moving mesh and a total Lagrangian one, on a case
# whose legacy arm eventually fails outright. Establishing which is closer to
# the truth needs a mesh study neither arm has had.
# ============================================================

echo "============================================================"
echo "ratCarotid regression test"
echo "============================================================"

failures=0

# Reference: the framework arm's own converged answer at the end time. The
# bounds are wide enough to survive a compiler or PETSc version change and
# narrow enough to catch the material or the formulation moving
MAG_D_MIN=3.2e-4
MAG_D_MAX=3.5e-4

run_framework() {
    local d="${REGRESSION_ROOT}/framework"

    rm -rf "${d}"; mkdir -p "${d}"
    for item in "${SCRIPT_DIR}"/*; do
        [[ "$(basename "${item}")" == "regressionTests" ]] && continue
        cp -a "${item}" "${d}/"
    done

    ( cd "${d}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true

    if solids4Foam::regressionCaseSkipped "${d}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: framework arm (the tutorial skipped here)"
        return 0
    fi

    if ! grep -q "Selecting mechanical constitutive law" \
        "${d}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the framework arm did not use the framework"
        return 1
    fi

    if ! grep -q "^End" "${d}/${SOLVER_LOGFILE}"; then
        echo "FAIL: the framework arm did not run to completion"
        tail -n 5 "${d}/${SOLVER_LOGFILE}" || true
        return 1
    fi

    if grep -qE "Nonlinear solve did not converge|SNES convergence error" \
        "${d}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the framework arm did not converge"
        return 1
    fi

    # The mixed formulation must have asked for a deviatoric implicit
    # stiffness. With the full one the bulk modulus - a thousand times the
    # shear modulus here - makes the Laplacian surrogate so stiff the linear
    # solve does not converge
    if ! grep -q "scalarDeviatoric" "${d}/${SOLVER_LOGFILE}"; then
        echo "FAIL: the mixed arm did not ask for a deviatoric stiffness"
        return 1
    fi

    local t m
    t=$(awk 'END {print $1}' "${d}/${DISP_FILE}" 2>/dev/null)
    m=$(awk 'END {print $5}' "${d}/${DISP_FILE}" 2>/dev/null)

    if [[ -z "${m}" ]]; then
        echo "FAIL: the framework arm produced no displacement history"
        return 1
    fi

    if ! awk "BEGIN {exit !((${t} - 1.0)^2 <= 1e-12)}"; then
        printf "FAIL: the framework arm stopped at t = %s, not the end time\n" \
            "${t}"
        return 1
    fi

    echo "PASS: the framework arm ran to completion and converged"

    if awk "BEGIN {exit !(${m} >= ${MAG_D_MIN} && ${m} <= ${MAG_D_MAX})}"; then
        printf "PASS: final inner-wall |D| = %.6g\n" "${m}"
    else
        printf "FAIL: final inner-wall |D| = %.6g (outside [%g, %g])\n" \
            "${m}" "${MAG_D_MIN}" "${MAG_D_MAX}"
        return 1
    fi

    # The law's own checks, which are what pins the constitutive port: an
    # honest isochoric split, and a fibre term that matches its closed form
    if command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
        local u="${REGRESSION_ROOT}/lawChecks"
        rm -rf "${u}"; mkdir -p "${u}"
        cp -a "${d}/constant" "${d}/system" "${u}/"
        rm -f "${u}/constant/solidProperties"
        cp -a "${d}/constant/solidProperties.framework" \
              "${u}/constant/solidProperties"

        # The closed-form fibre check needs the fibres along the stretch, so
        # the angle is zeroed here. It selects which directions the fibres
        # point in, not which code paths run
        sed -i.bak 's/\(fibreAngle.*\]\) *39.76;/\1 0.0;/' \
            "${u}/constant/mechanicalProperties"
        rm -f "${u}/constant/mechanicalProperties.bak"

        # Uniform directions, so the check needs no calcLocCoordinates run
        sed -i.bak 's|^        bulkModulus|        uniformLocalBasis yes;\n        Ec              (1 0 0);\n        Ea              (0 1 0);\n\n        bulkModulus|' \
            "${u}/constant/mechanicalProperties"
        rm -f "${u}/constant/mechanicalProperties.bak"

        if ! ( cd "${u}" && Test-mechanicalConstitutiveLaw > log.unit 2>&1 )
        then
            echo "FAIL: the law checks did not pass"
            grep -m3 "FAIL:" "${u}/log.unit" || true
            return 1
        fi

        if ! grep -q "fibre stress matches the closed form" "${u}/log.unit"
        then
            echo "FAIL: the fibre term was not checked against its closed form"
            return 1
        fi

        echo "PASS: $(grep -c 'PASS:' "${u}/log.unit") law checks, including" \
             "the fibre term against its closed form"
    fi

    return 0
}

if ! run_framework; then
    failures=$((failures + 1))
fi

echo
echo "NOTE: the legacy arm is not run. On foam-extend 4.1 it stalls and dies"
echo "      at t = 0.68, identically on the development branch, so there is"
echo "      no working reference to compare the ported law against."

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
