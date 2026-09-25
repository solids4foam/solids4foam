#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
elif command -v solids4FoamScripts.sh > /dev/null 2>&1; then
    source solids4FoamScripts.sh
fi

# Provide a fallback definition for older solids4FoamScripts.sh installs
# that pre-date the regressionCaseSkipped helper.
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
# idealisedVentricle regression test
# Runs the petsc and pressureDisplacement approaches with a
# shortened end-time and loosened tolerances, and checks that
# Max sigmaEq from the solver log is within sensible bounds.
# Approaches that are not supported in the current OpenFOAM
# flavour are reported as SKIP rather than failures.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

# The case is run to a small fraction of the full ramp time, so the resulting
# peak von Mises stress is well below the case maximum but still well above
# zero. It lands near 155 Pa - 154.725 at the time of writing - so the band
# below is wide enough to absorb the difference between forks and solvers
# while still being an oracle: [1, 1e7] admitted almost any number a run could
# produce, and so tested nothing on its own.
#
# The petsc arm is also required to stay within 1% of the removed legacy
# model's answer further down, which is the sharper of the two checks.
SIGMA_MIN=1.0e2
SIGMA_MAX=2.5e2

# The pressureDisplacement arm has a band of its own: it is a different mesh
# (the Fluent ventricle) and a fully incompressible material, and its peak von
# Mises stress at the same time is 719 Pa rather than 155 - so the band above
# never applied to it, which went unnoticed while it never ran to completion.
# 718.97 on foam-extend 4.1 with the loosened tolerances below and 718.78 with
# the tutorial's own, so +-5 % is loose on convergence and tight enough to
# catch the law or the formulation moving. There is no independent reference:
# at the first step the legacy law, which does not get further, differs by
# 4.5 % in this quantity (#466, #410)
SIGMA_PD_MIN=6.8e2
SIGMA_PD_MAX=7.6e2

# The petsc arm's Max sigmaEq on the removed legacy mechanicalModel, from the
# last commit that had it (mcl-stage8-coverage, c3a92b3d), on OpenFOAM.com
# v2512, the one fork the petsc arm runs on. The framework is not expected to
# match it exactly: its GuccioneElastic builds Q from the isochoric strain,
# where the legacy law built it from the full Green-Lagrange strain, so shape
# and volume are separated in one and coupled in the other. Both reduce to the
# published model in the incompressible limit it was written for. The 1% bound
# says the reformulation is the only thing between them; the framework reads
# 154.725 here, 0.17% below
LEGACY_PETSC_SIGMA=154.986
LEGACY_PETSC_REL_TOL=0.01

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

APPROACHES=(
    "petsc"
    "pressureDisplacement"
)

echo "============================================================"
echo "idealisedVentricle regression tests"
echo "Max sigmaEq in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "Approaches: ${APPROACHES[*]}"
echo "============================================================"
echo

# ------------------------------------------------------------
# Override helpers - shorten the case so the regression runs
# in a reasonable time without changing the published tutorial.
# ------------------------------------------------------------

shorten_controlDict() {
    local file="$1"
    if [[ -f "${file}" ]]; then
        sed -i.bak \
            -e 's/^\(\s*endTime\s*\).*/\1     0.02;/' \
            -e 's/^\(\s*deltaT\s*\).*/\1      0.01;/' \
            -e 's/^\(\s*writeControl\s*\).*/\1 timeStep;/' \
            -e 's/^\(\s*writeInterval\s*\).*/\1     1;/' \
            "${file}"
        rm -f "${file}.bak"
    fi
}

loosen_pressureDisplacement_tolerances() {
    local file="$1"
    if [[ -f "${file}" ]]; then
        sed -i.bak \
            -e 's/^\(\s*nCorrectors\s*\)[0-9]\+\s*;/\1            200;/' \
            -e 's/^\(\s*solutionTolerance\s*\)[0-9eE.+-]\+\s*;/\1      1e-04;/' \
            -e 's/^\(\s*alternativeTolerance\s*\)[0-9eE.+-]\+\s*;/\1   1e-03;/' \
            -e 's/^\(\s*materialTolerance\s*\)[0-9eE.+-]\+\s*;/\1      1e-03;/' \
            "${file}"
        rm -f "${file}.bak"
    fi
}

# There is no legacy answer for the pressureDisplacement arm to be compared
# with. The legacy GuccioneElastic pressureDisplacement mode did not get past
# the second step at any time step or tolerance tried, and did not before the
# framework existed either; the framework converges in a few hundred iterations
# a step (#466). The framework law takes its penalty from bulkModulus, which
# the tutorial sets to 1e15 so that the material is fully incompressible, as
# nu = 0.5 made it on the legacy path

prepare_case() {
    local case_dir="$1"
    local approach="$2"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        local base_item
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done

    # Shorten controlDict for the requested approach so the
    # regression run terminates after a couple of time steps.
    shorten_controlDict \
        "${case_dir}/caseOptions/${approach}/system/controlDict"

    if [[ "${approach}" == "pressureDisplacement" ]]; then
        loosen_pressureDisplacement_tolerances \
            "${case_dir}/caseOptions/${approach}/constant/solidProperties"
    fi
}

run_case() {
    local approach="$1"
    local case_dir="${REGRESSION_ROOT}/${approach}"

    prepare_case "${case_dir}" "${approach}"
    ( cd "${case_dir}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${case_dir}" && ./Allrun "${approach}" > "${ALLRUN_LOGFILE}" 2>&1 )
}

# ------------------------------------------------------------
# Parse args
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Run cases
# ------------------------------------------------------------

if [ "$CHECK_ONLY" = false ]; then
    rm -rf "${REGRESSION_ROOT}"
    mkdir -p "${REGRESSION_ROOT}"

    for approach in "${APPROACHES[@]}"; do
        run_case "${approach}"
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_sigma() {
    local case_dir="$1"

    grep "Max sigmaEq (von Mises stress)" "${case_dir}/${SOLVER_LOGFILE}" \
        2>/dev/null \
        | tail -n 1 \
        | awk '{print $NF}' \
        || true
}

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

for approach in "${APPROACHES[@]}"; do
    case_dir="${REGRESSION_ROOT}/${approach}"

    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "SKIP: ${approach} (not supported in this environment)"
        continue
    fi

    # The run has to have finished and converged. extract_max_sigma takes the
    # last value in the log, and a run that stopped early still leaves one -
    # so without this an arm that diverged at the second step would be
    # compared against the other arm's converged answer and could pass
    if ! grep -q "^End" "${case_dir}/${SOLVER_LOGFILE}"; then
        echo "FAIL: ${approach}: did not run to completion"
        failures=$((failures + 1))
        continue
    fi

    if grep -qE "Nonlinear solve did not converge|SNES convergence error" \
        "${case_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${approach}: did not converge"
        failures=$((failures + 1))
        continue
    fi

    if ! grep -q "Selecting mechanical constitutive law" \
        "${case_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${approach}: constructed no mechanical constitutive law"
        failures=$((failures + 1))
        continue
    fi

    sigma=$(extract_max_sigma "${case_dir}")

    if [[ -z "${sigma}" ]]; then
        echo "FAIL: ${approach}: could not extract Max sigmaEq"
        failures=$((failures + 1))
        continue
    fi

    sigma_min="${SIGMA_MIN}"
    sigma_max="${SIGMA_MAX}"
    if [[ "${approach}" == "pressureDisplacement" ]]; then
        sigma_min="${SIGMA_PD_MIN}"
        sigma_max="${SIGMA_PD_MAX}"
    fi

    if awk "BEGIN {exit !(${sigma} >= ${sigma_min} && ${sigma} <= ${sigma_max})}"
    then
        printf "PASS: %s: Max sigmaEq = %.6g\n" "${approach}" "${sigma}"
    else
        printf "FAIL: %s: Max sigmaEq = %.6g (outside [%g, %g])\n" \
            "${approach}" "${sigma}" "${sigma_min}" "${sigma_max}"
        failures=$((failures + 1))
    fi

    if [[ "${approach}" == "petsc" ]]; then
        petsc_sigma="${sigma}"
    fi
done

# The petsc arm against the removed legacy model, within the reformulation
if [[ -n "${petsc_sigma:-}" ]]; then
    if awk "BEGIN {exit !((${petsc_sigma} - ${LEGACY_PETSC_SIGMA})^2 \
        <= (${LEGACY_PETSC_REL_TOL}*${LEGACY_PETSC_SIGMA})^2)}"
    then
        printf "PASS: near the legacy model, differing by the reformulation (%.6g vs %.6g)\n" \
            "${petsc_sigma}" "${LEGACY_PETSC_SIGMA}"
    else
        printf "FAIL: differs from the legacy model by more than the reformulation explains (%.6g vs %.6g)\n" \
            "${petsc_sigma}" "${LEGACY_PETSC_SIGMA}"
        failures=$((failures + 1))
    fi
fi

# ------------------------------------------------------------
# The framework's check on GuccioneElastic's isochoric split
# ------------------------------------------------------------
# This case uses GuccioneElastic on its own, so its split can be checked
# directly: a law that declares it can separate its isochoric stress from its
# volumetric response is taken at its word by every mixed formulation, and a
# superposed dilation is the one thing that tells an honest split from dev()
# of a total. Where the law is wrapped in electroMechanicalLaw the check does
# not apply, because the active tension is not derived from a potential - so
# this is where it does apply
run_split_check() {
    # Run on any meshed arm that selected a mechanical constitutive law, so a
    # renamed or added arm remains covered, while an arm that stopped before
    # constructing its law is not selected merely because it appears first
    # (#466)
    local d
    for d in "${REGRESSION_ROOT}"/*; do
        [[ -d "${d}/constant/polyMesh" ]] || continue
        grep -q "Selecting mechanical constitutive law" \
            "${d}/${SOLVER_LOGFILE}" 2>/dev/null || continue

        if ! command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
            echo "SKIP: mechanicalConstitutiveLaw checks (not in PATH)"
            return 0
        fi

        if ! ( cd "${d}" && Test-mechanicalConstitutiveLaw > log.unit 2>&1 )
        then
            echo "FAIL: the law checks did not pass"
            grep -m2 "FAIL:" "${d}/log.unit" || true
            return 1
        fi

        if ! grep -q "isochoric stress ignores a superposed dilation" \
            "${d}/log.unit"
        then
            echo "FAIL: GuccioneElastic's isochoric split was not checked"
            return 1
        fi

        echo "PASS: GuccioneElastic's isochoric split is dilation invariant"
        return 0
    done

    echo "SKIP: mechanicalConstitutiveLaw checks (neither arm ran here)"
    return 0
}

if ! run_split_check; then
    failures=$((failures + 1))
fi

# ------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------

if [ "$CHECK_ONLY" = false ]; then
    for case_dir in "${REGRESSION_ROOT}"/*; do
        if [[ -d "${case_dir}" ]]; then
            ( cd "${case_dir}" && ./Allclean > /dev/null 2>&1 ) || true
        fi
    done
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
