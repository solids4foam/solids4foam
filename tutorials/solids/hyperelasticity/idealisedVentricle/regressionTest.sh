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

# The case is run to a small fraction of the full ramp time, so
# the resulting peak von Mises stress is well below the case
# maximum but still well above zero.
SIGMA_MIN=1.0e0
SIGMA_MAX=1.0e7

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

    sigma=$(extract_max_sigma "${case_dir}")

    if [[ -z "${sigma}" ]]; then
        echo "FAIL: ${approach}: could not extract Max sigmaEq"
        failures=$((failures + 1))
        continue
    fi

    if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${sigma} <= ${SIGMA_MAX})}"
    then
        printf "PASS: %s: Max sigmaEq = %.6g\n" "${approach}" "${sigma}"
    else
        printf "FAIL: %s: Max sigmaEq = %.6g (outside [%g, %g])\n" \
            "${approach}" "${sigma}" "${SIGMA_MIN}" "${SIGMA_MAX}"
        failures=$((failures + 1))
    fi
done

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
