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
# cantilever2d regression test
# Checks selected solution approaches against the analytical
# benchmark output.
# ============================================================

EPS_MIN=3.5e-4
EPS_MAX=5.0e-4
SIGMA_MIN=8.5e7
SIGMA_MAX=1.05e8
HIGH_ORDER_DISP_TOL=1e-10

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

APPROACHES=(
    petscSnes
    highOrder-movingLeastSquares
    highOrder-kExactLeastSquares
    highOrderJacobian
)

echo "============================================================"
echo "cantilever2d regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "High-order DDifference LInf < ${HIGH_ORDER_DISP_TOL}"
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

}

# Write a spatially varying initial stress into a case's 0 directory.
#
# Anything varying will do; the point is that a wrong cell-to-integration-point
# map changes the answer, which a uniform field would hide
write_sigma0_field() {
    local case_dir="$1"
    local nCells

    nCells=$(grep -h "nCells:" "${case_dir}/log.blockMesh" | tail -n 1 \
        | awk '{print $2}')

    if [[ -z "${nCells}" ]]; then
        echo "FAIL: could not read the cell count for the sigma0 field"
        return 1
    fi

    {
        echo "FoamFile"
        echo "{"
        echo "    version     2.0;"
        echo "    format      ascii;"
        echo "    class       volSymmTensorField;"
        echo "    location    \"0\";"
        echo "    object      sigma0;"
        echo "}"
        echo
        echo "dimensions      [1 -1 -2 0 0 0 0];"
        echo
        echo "internalField   nonuniform List<symmTensor>"
        echo "${nCells}"
        echo "("
        awk -v n="${nCells}" 'BEGIN {
            for (i = 0; i < n; i++)
            {
                f = 1e7*sin(0.01*i)
                printf "(%g %g %g %g %g %g)\n", \
                    f, 0.3*f, -0.2*f, 0.5*f, 0.1*f, -0.4*f
            }
        }'
        echo ")"
        echo ";"
        echo
        echo "boundaryField"
        echo "{"

        # A constraint patch - symmetry, empty, wedge - only takes the patch
        # field of its own kind, so these entries come from the mesh rather
        # than from a catch-all
        awk '
            /^\(/ { inList = 1; next }
            /^\)/ { inList = 0 }
            inList && /^[ \t]*[A-Za-z_][A-Za-z0-9_]*[ \t]*$/ {
                name = $1
                next
            }
            inList && /^[ \t]*type[ \t]/ {
                t = $2
                sub(";", "", t)
                if (name == "") next
                if (t == "patch" || t == "wall")
                {
                    t = "zeroGradient"
                }
                printf "    %s\n    {\n        type            %s;\n    }\n", \
                    name, t
                name = ""
            }
        ' "${case_dir}/constant/polyMesh/boundary"

        echo "}"
    } > "${case_dir}/0/sigma0"
}

# Run this case twice with an initial stress, once through the legacy
# mechanical law and once through the constitutive framework, and require the
# two to agree exactly.
#
# sigma0 is the first prescribed state the framework declares. It reaches the
# law by two routes - a uniform value in the law's dictionary, and a field the
# case supplies - and both are checked here because they are read by different
# code. A non-uniform field is used for the second so that a wrong
# cell-to-integration-point map changes the answer, which a uniform field
# would hide.
#
# This case is the host because it converges tightly enough for the two arms
# to agree to the last bit, which makes the check a statement about the model
# rather than about a solver tolerance
run_sigma0_comparison() {
    local mode="$1"
    local legacy_dir="${REGRESSION_ROOT}/sigma0Legacy-${mode}"
    local framework_dir="${REGRESSION_ROOT}/sigma0Framework-${mode}"
    local d

    for d in "${legacy_dir}" "${framework_dir}"; do
        rm -rf "${d}"
        mkdir -p "${d}"

        local item base_item
        for item in "${SCRIPT_DIR}"/*; do
            base_item=$(basename "${item}")
            if [[ "${base_item}" == "regressionTests" ]]; then
                continue
            fi
            cp -a "${item}" "${d}/"
        done

        if [[ "${mode}" == "dict" || "${mode}" == "both" ]]; then
            # A uniform initial stress given where the material is given
            sed -i \
                's|^\( *\)nu  *nu .*|&\n\1sigma0 sigma0 [1 -1 -2 0 0 0 0] (10e6 2e6 -3e6 15e6 0 -5e6);|' \
                "${d}/constant/mechanicalProperties"
        fi
    done

    # The two arms differ in this one entry and nothing else. It goes inside
    # the model's coeffs block, which is where the solid model looks for it;
    # at the top level it is read by nothing and silently ignored
    sed -i \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${framework_dir}/constant/solidProperties"

    for d in "${legacy_dir}" "${framework_dir}"; do
        (
            cd "${d}" || exit 1

            solids4Foam::convertCaseFormat . > log.convert 2>&1

            blockMesh > log.blockMesh 2>&1 || exit 1

            if [[ "${mode}" == "field" || "${mode}" == "both" ]]; then
                write_sigma0_field . || exit 1
            fi

            solids4Foam > log.solids4Foam 2>&1 || exit 1
        ) || { echo "FAIL: sigma0 ${mode} arm could not run"; return 1; }
    done

    # Each arm must have taken the path it was set up for, or the comparison
    # is between two copies of the same thing and proves nothing
    if grep -q "mechanicalConstitutiveLawManager" \
        "${legacy_dir}/log.solids4Foam"
    then
        echo "FAIL: the legacy sigma0 ${mode} arm used the framework"
        return 1
    fi

    if [[ "${mode}" == "dict" || "${mode}" == "both" ]]; then
        if ! grep -q "Uniform initial stress sigma0" \
            "${framework_dir}/log.solids4Foam"
        then
            echo "FAIL: the framework arm did not read sigma0 from the dict"
            return 1
        fi
    else
        if ! grep -q "Prescribed state 'sigma0'" \
            "${framework_dir}/log.solids4Foam"
        then
            echo "FAIL: the framework arm did not read the sigma0 field"
            return 1
        fi
    fi

    local t
    t=$(foamListTimes -case "${legacy_dir}" -latestTime 2>/dev/null | tail -n 1)

    if [[ -z "${t}" ]]; then
        echo "FAIL: sigma0 ${mode} arms produced no result"
        return 1
    fi

    if [[ ! -f "${legacy_dir}/${t}/D" || ! -f "${framework_dir}/${t}/D" ]]; then
        echo "FAIL: sigma0 ${mode} arms produced no D field"
        return 1
    fi

    # sigma0 must actually have changed the answer, or agreement is vacuous
    if diff -q "${legacy_dir}/${t}/D" "${SIGMA0_BASELINE_D}" > /dev/null 2>&1
    then
        echo "FAIL: sigma0 ${mode} left the solution unchanged"
        return 1
    fi

    if ! diff -q "${legacy_dir}/${t}/D" "${framework_dir}/${t}/D" > /dev/null
    then
        echo "FAIL: sigma0 ${mode} legacy and framework differ"
        return 1
    fi

    # Legacy reads a sigma0 field and then assigns any dictionary sigma0 over
    # the whole of it, so when a case carries both, the dictionary is what
    # takes effect. Checking that here is the only way to see that the
    # framework resolves the two the same way round
    if [[ "${mode}" == "dict" ]]; then
        SIGMA0_DICT_D="${legacy_dir}/${t}/D"
    elif [[ "${mode}" == "both" ]]; then
        if ! diff -q "${legacy_dir}/${t}/D" "${SIGMA0_DICT_D}" > /dev/null
        then
            echo "FAIL: sigma0 both: the field was not overridden by the dict"
            return 1
        fi
    fi

    echo "PASS: sigma0 ${mode} legacy and framework agree exactly"
    return 0
}

# A prescribed field has to survive a restart. It is written into 0 and the
# case then restarts from a later time whose directory does not contain it, so
# this is the check that the field is looked for where it actually lives.
#
# It compares the framework against itself, restarted against continuous,
# rather than against legacy. Legacy looks for sigma0 only beside the fields
# it restarts from and so loses it, keeping it only when an earlier run had
# already written it forward. That is a hole rather than a behaviour worth
# reproducing, so the framework diverges here deliberately
run_sigma0_restart_check() {
    local root="${REGRESSION_ROOT}"
    local d

    for d in restartContinuous restartRestarted restartNone; do
        rm -rf "${root}/${d}"
        mkdir -p "${root}/${d}"

        local item base_item
        for item in "${SCRIPT_DIR}"/*; do
            base_item=$(basename "${item}")
            if [[ "${base_item}" == "regressionTests" ]]; then
                continue
            fi
            cp -a "${item}" "${root}/${d}/"
        done

        sed -i \
            's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
            "${root}/${d}/constant/solidProperties"
    done

    for d in restartContinuous restartRestarted restartNone; do
        (
            cd "${root}/${d}" || exit 1
            solids4Foam::convertCaseFormat . > log.convert 2>&1
            blockMesh > log.blockMesh 2>&1 || exit 1
        ) || { echo "FAIL: sigma0 restart check could not mesh ${d}"; return 1; }
    done

    # The two that carry an initial stress get it in 0, and only in 0
    for d in restartContinuous restartRestarted; do
        write_sigma0_field "${root}/${d}" || return 1
    done

    (
        cd "${root}/restartContinuous" || exit 1
        sed -i 's|^endTime .*|endTime         2;|' system/controlDict
        solids4Foam > log.solids4Foam 2>&1
    ) || { echo "FAIL: sigma0 restart check could not run continuous"; return 1; }

    (
        cd "${root}/restartNone" || exit 1
        sed -i 's|^endTime .*|endTime         2;|' system/controlDict
        solids4Foam > log.solids4Foam 2>&1
    ) || { echo "FAIL: sigma0 restart check could not run the control"; return 1; }

    (
        cd "${root}/restartRestarted" || exit 1
        solids4Foam > log.first 2>&1 || exit 1
        sed -i \
            's|^endTime .*|endTime         2;|;s|^startFrom .*|startFrom       latestTime;|' \
            system/controlDict
        solids4Foam > log.solids4Foam 2>&1
    ) || { echo "FAIL: sigma0 restart check could not restart"; return 1; }

    # The restart has to be a real one: the time it resumes from must not
    # itself contain sigma0, or the fallback is never exercised
    if [[ -f "${root}/restartRestarted/1/sigma0" ]]; then
        echo "FAIL: sigma0 restart check did not test a restart"
        return 1
    fi

    local cont="${root}/restartContinuous/2/D"
    local rest="${root}/restartRestarted/2/D"
    local none="${root}/restartNone/2/D"

    if [[ ! -f "${cont}" || ! -f "${rest}" || ! -f "${none}" ]]; then
        echo "FAIL: sigma0 restart check produced no D field"
        return 1
    fi

    if diff -q "${cont}" "${none}" > /dev/null; then
        echo "FAIL: sigma0 restart check: sigma0 changed nothing"
        return 1
    fi

    if diff -q "${cont}" "${rest}" > /dev/null; then
        echo "PASS: sigma0 survives a restart"
        return 0
    fi

    echo "FAIL: sigma0 was lost on restart"
    return 1
}

# A run with no initial stress, to show that the ones with it are different
make_sigma0_baseline() {
    local d="${REGRESSION_ROOT}/sigma0Baseline"
    local item base_item

    rm -rf "${d}"
    mkdir -p "${d}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${d}/"
    done

    (
        cd "${d}" || exit 1
        solids4Foam::convertCaseFormat . > log.convert 2>&1
        blockMesh > log.blockMesh 2>&1 || exit 1
        solids4Foam > log.solids4Foam 2>&1 || exit 1
    ) || return 1

    local t
    t=$(foamListTimes -case "${d}" -latestTime 2>/dev/null | tail -n 1)
    SIGMA0_BASELINE_D="${d}/${t}/D"
}

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
    prepare_case
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

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

extract_disp_linf() {
    grep -A2 "Writing DDifference field" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | grep "Norms:" -A1 \
        | tail -n 1 \
        | awk '{print $3}' \
        || true
}

select_run_approach() {
    local requested="$1"

    case "${requested}" in
        highOrder-movingLeastSquares|highOrder-kExactLeastSquares)
            local least_squares_type="${requested#highOrder-}"
            sed -E -i.bak \
                "s/^([[:space:]]*)type[[:space:]]+(movingLeastSquares|kExactLeastSquares);/\\1type ${least_squares_type};/" \
                "${CASE_DIR}/constant/solidProperties.highOrder"
            rm -f "${CASE_DIR}/constant/solidProperties.highOrder.bak"
            RUN_APPROACH=highOrder
            ;;
        *)
            RUN_APPROACH="${requested}"
            ;;
    esac
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

check_high_order_errors() {
    local displacement
    local failures=0

    displacement=$(extract_disp_linf)

    if [[ -z "${displacement}" ]]; then
        echo "FAIL: Could not extract high-order DDifference LInf"
        return 1
    fi

    if awk "BEGIN {exit !(${displacement} < ${HIGH_ORDER_DISP_TOL})}"; then
        printf "PASS: High-order DDifference LInf = %.6g\n" "${displacement}"
    else
        printf "FAIL: High-order DDifference LInf = %.6g\n" "${displacement}"
        failures=$((failures + 1))
    fi

    return "${failures}"
}

failures=0

if [ "$CHECK_ONLY" = false ]; then
    for approach in "${APPROACHES[@]}"; do
        echo
        echo "------------------------------------------------------------"
        echo "Testing approach: ${approach}"
        echo "------------------------------------------------------------"

        select_run_approach "${approach}"
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
        ( cd "${CASE_DIR}" && ./Allrun "${RUN_APPROACH}" > "${ALLRUN_LOGFILE}" 2>&1 )

        if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
            echo "Skipping ${approach} because it is unavailable in this environment"
            continue
        fi

        if ! check_solver_extrema "${approach}"; then
            failures=$((failures + 1))
        fi

        if [[ "${RUN_APPROACH}" == highOrder* ]] && ! check_high_order_errors; then
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

    if grep -q "highOrderResidual true" "${CASE_DIR}/constant/solidProperties" \
        && ! check_high_order_errors; then
        failures=$((failures + 1))
    fi
fi

if [ "$CHECK_ONLY" = false ]; then
    if make_sigma0_baseline; then
        for sigma0_mode in dict field both; do
            if ! run_sigma0_comparison "${sigma0_mode}"; then
                failures=$((failures + 1))
            fi
        done

        if ! run_sigma0_restart_check; then
            failures=$((failures + 1))
        fi
    else
        echo "FAIL: could not run the sigma0 baseline"
        failures=$((failures + 1))
    fi
fi

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
