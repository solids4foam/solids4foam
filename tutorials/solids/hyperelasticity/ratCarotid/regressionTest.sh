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
# Two formulations, both on the mechanicalConstitutiveLaw framework, and only
# one of them currently completes.
#
# coupledPressureDisplacementSolid (./Allrun pressureDisplacement) is
# foam-extend only and does not reach the end time: it stalls at a relative
# residual near 1 and stops at t = 0.3, when the deformation gradient inverts.
# That is longstanding rather than a regression - with the former legacy
# HolzapfelGasserOgdenElastic it stalled the same way and stopped at t = 0.68.
# It is not run here.
#
# nonLinearGeometryTotalLagrangianTotalDisplacement with solvePressure, the
# default, runs to completion on any fork, which is what this checks.
#
# It requires pressure stabilisation. With momentum stabilisation disabled it
# reaches t = 0.12 before the nonlinear solve stalls. A representative
# momentum scale of 1 completes the case; the previous value of 100 was
# unnecessarily large. These terms affect the discrete equations, and the
# legacy law never completed the case, so there is no legacy answer to compare
# with: this test pins the framework result and verifies the constitutive-law
# checks.
# ============================================================

echo "============================================================"
echo "ratCarotid regression test"
echo "============================================================"

failures=0

# Reference: the case's own converged answer at the end time. The
# bounds are wide enough to survive a compiler or PETSc version change and
# narrow enough to catch the material or the formulation moving
MAG_D_MIN=3.2e-4
MAG_D_MAX=3.5e-4

run_case() {
    local d="${REGRESSION_ROOT}/main"

    rm -rf "${d}"; mkdir -p "${d}"
    for item in "${SCRIPT_DIR}"/*; do
        [[ "$(basename "${item}")" == "regressionTests" ]] && continue
        cp -a "${item}" "${d}/"
    done

    ( cd "${d}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true

    if solids4Foam::regressionCaseSkipped "${d}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: the tutorial skipped here"
        return 0
    fi

    if ! grep -q "Selecting mechanical constitutive law" \
        "${d}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the case constructed no mechanical constitutive law"
        return 1
    fi

    if ! grep -q "^End" "${d}/${SOLVER_LOGFILE}"; then
        echo "FAIL: the case did not run to completion"
        tail -n 5 "${d}/${SOLVER_LOGFILE}" || true
        return 1
    fi

    if grep -qE "Nonlinear solve did not converge|SNES convergence error" \
        "${d}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the case did not converge"
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
        echo "FAIL: the case produced no displacement history"
        return 1
    fi

    if ! awk "BEGIN {exit !((${t} - 1.0)^2 <= 1e-12)}"; then
        printf "FAIL: the case stopped at t = %s, not the end time\n" \
            "${t}"
        return 1
    fi

    echo "PASS: the case ran to completion and converged"

    if awk "BEGIN {exit !(${m} >= ${MAG_D_MIN} && ${m} <= ${MAG_D_MAX})}"; then
        printf "PASS: final inner-wall |D| = %.6g\n" "${m}"
    else
        printf "FAIL: final inner-wall |D| = %.6g (outside [%g, %g])\n" \
            "${m}" "${MAG_D_MIN}" "${MAG_D_MAX}"
        return 1
    fi

    # The law's own checks, which are what pins the constitutive port: an
    # honest isochoric split, and a fibre term that matches its closed form
    local testApp=0
    solids4Foam::requireTestApp Test-mechanicalConstitutiveLaw || testApp=$?

    # A failure in CI, where the application is always built
    if (( testApp == 2 )); then
        return 1
    fi

    if (( testApp == 0 )); then
        local u="${REGRESSION_ROOT}/lawChecks"
        rm -rf "${u}"; mkdir -p "${u}"
        cp -a "${d}/constant" "${d}/system" "${u}/"
        rm -f "${u}/constant/solidProperties"
        cp -a "${d}/constant/solidProperties.totalLagrangian" \
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

if ! run_case; then
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
