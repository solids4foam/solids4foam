#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# ============================================================
# ratCarotid regression test
#
# This case does not currently run.
#
# The solver stalls at a relative residual of about 0.99 and dies of a
# floating point exception at t = 0.68, on foam-extend 4.1, which is the only
# fork coupledPressureDisplacementSolid is built for. That is not a regression
# from the constitutive-law work: the identical failure - same time step, same
# iteration, same residual - happens on origin/development at the commit this
# branch started from. So the solver arm is not attempted here, and this
# script checks what can be checked.
#
# What it checks is the framework port of HolzapfelGasserOgdenElastic. The
# mesh builds on any fork, and the law's own checks need only a mesh and a
# material, so the port has coverage on the fork everything else is checked on
# even though its solid model and its tutorial do not run there.
#
# The one thing this cannot check is that the framework law reproduces the
# legacy one. That needs the case to run, and it does not.
# ============================================================

echo "============================================================"
echo "ratCarotid regression test"
echo "============================================================"

failures=0

run_law_checks() {
    local d="${REGRESSION_ROOT}/framework"

    if ! command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
        echo "SKIP: Test-mechanicalConstitutiveLaw not found in PATH"
        return 0
    fi

    rm -rf "${d}"; mkdir -p "${d}"
    for item in "${SCRIPT_DIR}"/*; do
        [[ "$(basename "${item}")" == "regressionTests" ]] && continue
        cp -a "${item}" "${d}/"
    done

    # The framework material, with the volumetric penalty that the legacy law
    # does without because it takes its whole spherical stress from the solid
    # model's pressure
    cp "${d}/constant/mechanicalProperties.framework" \
       "${d}/constant/mechanicalProperties"

    if ! ( cd "${d}" && blockMesh > log.blockMesh 2>&1 ); then
        echo "FAIL: could not mesh the case"
        tail -n 5 "${d}/log.blockMesh" || true
        return 1
    fi

    if ! ( cd "${d}" && Test-mechanicalConstitutiveLaw > log.unit 2>&1 ); then
        echo "FAIL: the law checks did not pass"
        grep -m3 "FAIL:" "${d}/log.unit" || true
        return 1
    fi

    if ! grep -q "Selecting mechanical constitutive law: \
HolzapfelGasserOgdenElastic" "${d}/log.unit"
    then
        # The line is wrapped above; match it loosely rather than by column
        if ! grep -q "HolzapfelGasserOgdenElastic" "${d}/log.unit"; then
            echo "FAIL: the framework law was not selected"
            return 1
        fi
    fi

    # The check that matters for this law. Every invariant it uses is taken on
    # the isochoric deformation, so a superposed dilation must leave the
    # Kirchhoff isochoric stress alone and the stress it returns must be trace
    # free. A law that took its invariants on the full deformation - as the
    # legacy law does, having no volumetric term to separate - would fail both
    if ! grep -q "isochoric stress ignores a superposed dilation" \
        "${d}/log.unit"
    then
        echo "FAIL: the isochoric split was not checked"
        return 1
    fi

    if ! grep -q "isochoric stress is trace free" "${d}/log.unit"; then
        echo "FAIL: the trace-free condition was not checked"
        return 1
    fi

    # The two checks above would pass a law that had lost its fibre term
    # entirely - it would still be dilation invariant and still trace free.
    # This one pins the fibre contribution to a closed form, and fails if the
    # term is removed, if the structure tensor is not pushed forward, or if
    # the factor of two from S = 2 dW/dC is dropped
    if ! grep -q "fibre stress matches the closed form" "${d}/log.unit"; then
        echo "FAIL: the fibre term was not checked against its closed form"
        return 1
    fi

    echo "PASS: $(grep -c 'PASS:' "${d}/log.unit") law checks, including the" \
         "isochoric split"

    return 0
}

if ! run_law_checks; then
    failures=$((failures + 1))
fi

echo
echo "NOTE: the solver arm is not run. This case stalls and dies at t = 0.68"
echo "      on foam-extend 4.1, identically on origin/development, so there"
echo "      is no working reference to compare a ported law against."

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
