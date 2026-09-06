#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"
SOLIDS4FOAM_ROOT_ABS=$(cd "${SCRIPT_DIR}/../../../../" && pwd)

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# sphericalCavity regression test
# Checks selected solution approaches against the expected solution bounds
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

EPS_MIN=7.0e-6
EPS_MAX=9.0e-6
SIGMA_MIN=1.8e6
SIGMA_MAX=2.0e6

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
PARALLEL_N_PROCS=2

APPROACHES=(
    segregated
    petscSnes
    highOrder-movingLeastSquares
    highOrder-kExactLeastSquares
    highOrder-movingLeastSquares-parallel
    highOrder-kExactLeastSquares-parallel
)

echo "============================================================"
echo "sphericalCavity regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
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

    # The regression copy lives deeper than the source tutorial, so the
    # relative SOLIDS4FOAM_ROOT in this local library build no longer points to
    # the repository root.
    sed -i.bak \
        "s|^SOLIDS4FOAM_ROOT := .*|SOLIDS4FOAM_ROOT := ${SOLIDS4FOAM_ROOT_ABS}|" \
        "${CASE_DIR}/src/Make/options"

    sed -E -i.bak \
        "s/^[[:space:]]*numberOfSubdomains[[:space:]]+[0-9]+;/numberOfSubdomains ${PARALLEL_N_PROCS};/" \
        "${CASE_DIR}/system/decomposeParDict"
    rm -f "${CASE_DIR}/system/decomposeParDict.bak"

}

# ------------------------------------------------------------
# Clean & run case
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

if [ "$CHECK_ONLY" = false ]; then
    if ! command -v gmsh > /dev/null 2>&1; then
        echo "Skipping regression checks because Gmsh is not installed"
        exit 0
    fi

    prepare_case
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null \
        | tail -n 1 \
        | awk '{print $NF}' || true
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null \
        | tail -n 1 \
        | awk '{print $NF}' || true
}

select_run_approach() {
    local requested="$1"

    case "${requested}" in
        highOrder-movingLeastSquares|highOrder-kExactLeastSquares)
            local least_squares_type="${requested#highOrder-}"
            sed -E -i.bak \
                -e "s/^([[:space:]]*)type[[:space:]]+(movingLeastSquares|kExactLeastSquares);/\\1type ${least_squares_type};/" \
                -e "s/^([[:space:]]*)highOrderJacobian[[:space:]]+(true|false);/\\1highOrderJacobian false;/" \
                "${CASE_DIR}/constant/solidProperties.highOrder"
            rm -f "${CASE_DIR}/constant/solidProperties.highOrder.bak"
            RUN_APPROACH=highOrder
            RUN_MODE=serial
            ;;
        highOrder-movingLeastSquares-parallel|highOrder-kExactLeastSquares-parallel)
            local least_squares_type="${requested#highOrder-}"
            least_squares_type="${least_squares_type%-parallel}"
            sed -E -i.bak \
                -e "s/^([[:space:]]*)type[[:space:]]+(movingLeastSquares|kExactLeastSquares);/\\1type ${least_squares_type};/" \
                -e "s/^([[:space:]]*)highOrderJacobian[[:space:]]+(true|false);/\\1highOrderJacobian false;/" \
                "${CASE_DIR}/constant/solidProperties.highOrder"
            rm -f "${CASE_DIR}/constant/solidProperties.highOrder.bak"
            RUN_APPROACH=highOrder
            RUN_MODE=parallel
            ;;
        *)
            RUN_APPROACH="${requested}"
            RUN_MODE=serial
            ;;
    esac
}


link_case_files_for_suffix() {
    local suffix="$1"

    while IFS= read -r -d '' file; do
        rm -f "${file%.*}"
        ln -nsf "$(basename "${file}")" "${file%.*}"
    done < <(
        find "${CASE_DIR}/constant" "${CASE_DIR}/system" \
            -type f -name "*.${suffix}" -print0
    )
}


run_parallel_case() {
    solids4Foam::caseDoesNotRunWithOpenFOAMOrg
    solids4Foam::requirePetscOrExitSilently

    link_case_files_for_suffix highOrder

    (
        cd "${CASE_DIR}"

        echo "Compiling libraries..."
        (cd src && bash ./Allwmake -s)

        mv 0 0.tmp
        solids4Foam::runApplication gmsh -3 -format msh2 sphericalCavity.geo
        solids4Foam::runApplication gmshToFoam sphericalCavity.msh
        solids4Foam::runApplication changeDictionary
        mv 0.tmp 0

        solids4Foam::runApplication decomposePar
        solids4Foam::runParallel solids4Foam
        solids4Foam::runApplication reconstructPar
    )
}

run_parallel_high_order_grad_test() {
    local least_squares_type="$1"
    local suffix="${least_squares_type}.parallel"
    local log_file="${CASE_DIR}/log.Test-highOrderGrad.${suffix}"

    if ! command -v Test-highOrderGrad >/dev/null 2>&1; then
        echo "SKIP: Test-highOrderGrad is not available"
        return 0
    fi

    sed -E -i.bak \
        "s/^([[:space:]]*)type[[:space:]]+(movingLeastSquares|kExactLeastSquares);/\1type ${least_squares_type};/" \
        "${CASE_DIR}/constant/solidProperties.highOrder"
    rm -f "${CASE_DIR}/constant/solidProperties.highOrder.bak"

    if ! grep -qE \
        "^[[:space:]]*type[[:space:]]+${least_squares_type};" \
        "${CASE_DIR}/constant/solidProperties.highOrder"
    then
        echo "FAIL: Could not select ${least_squares_type}"
        return 1
    fi

    if ! (
        cd "${CASE_DIR}"
        solids4Foam::runParallel \
            -n "${PARALLEL_N_PROCS}" \
            -o -s "${suffix}" Test-highOrderGrad >/dev/null 2>&1
    ); then
        echo "FAIL: Test-highOrderGrad (${least_squares_type}, parallel)"
        return 1
    fi

    if grep -q "Overall result: PASSED" "${log_file}"
    then
        echo "PASS: Test-highOrderGrad (${least_squares_type}, parallel)"
        return 0
    fi

    echo "FAIL: Test-highOrderGrad (${least_squares_type}, parallel)"
    return 1
}


check_solver_extrema() {
    local approach="$1"
    local epsilon
    local sigma
    local failures=0

    epsilon=$(extract_max_epsilon)
    sigma=$(extract_max_sigma)

    if [[ -z "${epsilon}" || -z "${sigma}" ]]; then
        echo "FAIL: Could not extract one or more regression quantities for ${approach}"
        return 1
    fi

    if awk "BEGIN {exit !(${epsilon} >= ${EPS_MIN} && ${epsilon} <= ${EPS_MAX})}"; then
        printf "PASS: Max epsilonEq = %.6g\n" "${epsilon}"
    else
        printf "FAIL: Max epsilonEq = %.6g\n" "${epsilon}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${sigma} <= ${SIGMA_MAX})}"; then
        printf "PASS: Max sigmaEq = %.6g\n" "${sigma}"
    else
        printf "FAIL: Max sigmaEq = %.6g\n" "${sigma}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if [ "$CHECK_ONLY" = false ]; then
    for approach in "${APPROACHES[@]}"; do
        echo
        echo "------------------------------------------------------------"
        echo "Testing approach: ${approach}"
        echo "------------------------------------------------------------"

        select_run_approach "${approach}"
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

        if [[ "${RUN_MODE}" == "parallel" ]]; then
            run_parallel_case \
                > "${CASE_DIR}/${ALLRUN_LOGFILE}" 2>&1
        else
            (
                cd "${CASE_DIR}"
                ./Allrun "${RUN_APPROACH}" \
                    > "${ALLRUN_LOGFILE}" 2>&1
            )
        fi

        if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
            echo "Skipping ${approach} because it is unavailable in this environment"
            continue
        fi

        if [[ "${approach}" == highOrder-*-parallel ]]; then
            least_squares_type="${approach#highOrder-}"
            least_squares_type="${least_squares_type%-parallel}"

            if ! run_parallel_high_order_grad_test "${least_squares_type}"; then
                failures=$((failures + 1))
            fi
        fi

        if ! check_solver_extrema "${approach}"; then
            failures=$((failures + 1))
        fi
    done
else
    if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    if ! check_solver_extrema "check-only"; then
        failures=$((failures + 1))
    fi
fi

# Clean case again
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
