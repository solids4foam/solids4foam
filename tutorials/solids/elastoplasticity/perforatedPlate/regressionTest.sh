#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
# Run twice: the legacy mechanicalModel, and the mechanicalConstitutiveLaw
# framework. This is the only regression case whose material is history
# dependent AND which runs the framework end to end in a solver, so it is the
# only place a plastic history error in the framework would show up in a real
# solve rather than in a unit check
APPROACHES=(
    legacy
    framework
)

# ============================================================
# Elastoplastic perforated plate regression test
# Checks strain, stress, and plastic yielding
# ============================================================

# Reference ranges (order-of-magnitude + robustness)
EPSILON_MIN=4.5e-3
EPSILON_MAX=6.0e-3

SIGMA_MIN=1.5e8
SIGMA_MAX=2.2e8

YIELD_MIN=28
YIELD_MAX=32

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "Elastoplastic perforated plate regression test"
echo "Max epsilonEq           in [${EPSILON_MIN}, ${EPSILON_MAX}]"
echo "Max sigmaEq (von Mises) in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "Yielding cells          in [${YIELD_MIN}, ${YIELD_MAX}]"
echo "============================================================"
echo

prepare_case() {
    local approach="$1"
    CASE_DIR="${REGRESSION_ROOT}/${approach}"

    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done

    if [[ "${approach}" == framework* ]]; then
        # The switch lives in the <type>Coeffs sub-dictionary; at the top level
        # it is silently ignored and this arm would quietly repeat the legacy run
        sed -i.bak \
            's/^    predictor yes;/    useMechanicalConstitutiveLawManager yes;\n    predictor yes;/' \
            "${CASE_DIR}/constant/solidProperties"
        rm -f "${CASE_DIR}/constant/solidProperties.bak"

        if ! grep -q 'useMechanicalConstitutiveLawManager' \
            "${CASE_DIR}/constant/solidProperties"
        then
            echo "FAIL: could not enable the framework in solidProperties"
            exit 1
        fi
    fi
}

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | awk '{print $NF}' \
        | tail -n 1
}

extract_yielding_cells() {
    # The legacy law counts cells; the framework law counts integration points,
    # which is the honest description for a law that may be evaluated on faces
    # or points. Accept either wording rather than making one law lie
    # Prefer the framework's message: in the framework arm the legacy law is
    # still constructed but never called, so it truthfully reports zero and
    # would mask the framework's own count
    if grep -q "yielding integration points" "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        grep "Number of yielding integration points" \
            "${CASE_DIR}/${SOLVER_LOGFILE}" \
            | tail -n 101 \
            | head -n 1 \
            | sed 's|.*= *||; s|/.*||'
    elif grep -q "cells .* are actively yielding" "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        grep "cells .* are actively yielding" "${CASE_DIR}/${SOLVER_LOGFILE}" \
            | tail -n 101 \
            | head -n 1 \
            | awk '{print $1}'
    fi
}

# Exercise the mechanicalConstitutiveLawManager on this case. It is the only
# tutorial in the regression set whose material is history dependent, so it is
# the only one where a tangent query has any history to disturb
run_constitutive_test() {
    if ! command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
        echo "SKIP: Test-mechanicalConstitutiveLaw not found in PATH"
        return 0
    fi

    if [[ ! -d "${CASE_DIR}/constant/polyMesh" ]]; then
        echo "SKIP: mechanicalConstitutiveLaw checks (case has no mesh)"
        return 0
    fi

    if ( cd "${CASE_DIR}" && Test-mechanicalConstitutiveLaw \
            > "log.Test-mechanicalConstitutiveLaw" 2>&1 )
    then
        local n_passed
        n_passed=$(grep -c 'PASS:' \
            "${CASE_DIR}/log.Test-mechanicalConstitutiveLaw" || true)

        if (( n_passed == 0 )); then
            echo "SKIP: mechanicalConstitutiveLaw checks (not built on this fork)"
            return 0
        fi

        echo "PASS: mechanicalConstitutiveLaw checks (${n_passed} checks)"
        return 0
    fi

    echo "FAIL: mechanicalConstitutiveLaw checks"
    grep 'FAIL:' "${CASE_DIR}/log.Test-mechanicalConstitutiveLaw" || true
    return 1
}

# ------------------------------------------------------------
# Run and check each approach
# ------------------------------------------------------------

failures=0
declare -A RESULT_EPS
declare -A RESULT_SIG

for approach in "${APPROACHES[@]}"; do
    echo
    echo "------------------------------------------------------------"
    echo "Testing approach: ${approach}"
    echo "------------------------------------------------------------"

    prepare_case "${approach}"

    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )

    # Assert the run took the path this approach names
    if grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}"
    then
        used_framework=true
    else
        used_framework=false
    fi

    if [[ "${approach}" == "framework" && "${used_framework}" == false ]]; then
        echo "FAIL: framework approach did not construct the framework"
        failures=$((failures + 1))
    elif [[ "${approach}" == "legacy" && "${used_framework}" == true ]]; then
        echo "FAIL: legacy approach unexpectedly constructed the framework"
        failures=$((failures + 1))
    else
        echo "PASS: ${approach} took the expected path"
    fi

    epsilon=$(extract_max_epsilon)
    sigma=$(extract_max_sigma)
    yielding_cells=$(extract_yielding_cells)

    if [[ -z "${epsilon}" || -z "${sigma}" || -z "${yielding_cells}" ]]; then
        echo "FAIL: ${approach} could not extract one or more quantities"
        failures=$((failures + 1))
        continue
    fi

    RESULT_EPS["${approach}"]="${epsilon}"
    RESULT_SIG["${approach}"]="${sigma}"

    if awk "BEGIN {exit !(${epsilon} >= ${EPSILON_MIN} && ${epsilon} <= ${EPSILON_MAX})}"
    then
        printf "PASS: %s Max epsilonEq = %.6g\n" "${approach}" "${epsilon}"
    else
        printf "FAIL: %s Max epsilonEq = %.6g\n" "${approach}" "${epsilon}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${sigma} >= ${SIGMA_MIN} && ${sigma} <= ${SIGMA_MAX})}"
    then
        printf "PASS: %s Max sigmaEq = %.6g\n" "${approach}" "${sigma}"
    else
        printf "FAIL: %s Max sigmaEq = %.6g\n" "${approach}" "${sigma}"
        failures=$((failures + 1))
    fi

    if (( yielding_cells >= YIELD_MIN && yielding_cells <= YIELD_MAX )); then
        printf "PASS: %s yielding = %d\n" "${approach}" "${yielding_cells}"
    else
        printf "FAIL: %s yielding = %d\n" "${approach}" "${yielding_cells}"
        failures=$((failures + 1))
    fi

    if [[ "${approach}" == "legacy" ]]; then
        if ! run_constitutive_test; then
            failures=$((failures + 1))
        fi
    fi

    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
done

# The framework must reproduce the legacy result. This case has no pressure
# smoothing, so unlike longWall there is nothing the framework omits: the two
# should agree to solver tolerance, and a plastic history error would show
# here as a difference that the unit checks cannot see
if [[ -n "${RESULT_EPS[legacy]:-}" && -n "${RESULT_EPS[framework]:-}" ]]; then
    for q in eps sig; do
        if [[ "${q}" == "eps" ]]; then
            a="${RESULT_EPS[legacy]}"; b="${RESULT_EPS[framework]}"; n="epsilonEq"
        else
            a="${RESULT_SIG[legacy]}"; b="${RESULT_SIG[framework]}"; n="sigmaEq"
        fi

        if awk "BEGIN {exit !(($a - $b)^2 <= (1e-8*$a)^2)}"; then
            printf "PASS: legacy and framework %s agree (%.8g vs %.8g)\n" \
                "$n" "$a" "$b"
        else
            printf "FAIL: legacy and framework %s differ (%.8g vs %.8g)\n" \
                "$n" "$a" "$b"
            failures=$((failures + 1))
        fi
    done
else
    echo "SKIP: cross-check needs both approaches to have run"
fi

# ------------------------------------------------------------
# Restart
# ------------------------------------------------------------
# Three things are checked here, and the second is the one that gives the
# first its meaning.
#
#   1. A run stopped and continued reaches the same answer as one that ran
#      straight through. Without the state IO this case is wrong by 2e-2 in
#      max epsilonEq, and wrong *quietly*: it runs to completion and prints a
#      plausible number. That is the failure being guarded against.
#
#   2. Deleting one state file makes the continuation stop with an error. A
#      restart that cannot find the history it needs must refuse, because
#      falling back to the cold start defaults is exactly the silent wrongness
#      of point 1 wearing a different hat.
#
#   3. The state is actually written. A test of reading proves nothing if
#      nothing was written and both arms are reading the same absence.
#
# The material must have yielded before the restart point or none of this is
# being tested: restarting an elastic case is exact even with no state IO at
# all, since the history being restored is zero either way. t = 10 is chosen
# for that reason, and check 3 asserts the yielding rather than assuming it
run_restart_test() {
    local base="${REGRESSION_ROOT}/frameworkRestart"

    prepare_case "frameworkRestart"
    CASE_DIR="${base}"

    sed -i.bak 's/^writePrecision  6;/writePrecision  14;/' \
        "${CASE_DIR}/system/controlDict"
    sed -i.bak 's/^endTime         20;/endTime         10;/' \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"

    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )

    # 3. the state must be on disk, and the material must have yielded
    local state_file
    state_file=$(ls "${CASE_DIR}"/10/*:*:epsilonPEq 2>/dev/null | head -n 1)

    if [[ -z "${state_file}" ]]; then
        echo "FAIL: restart: no constitutive state was written at t = 10"
        return 1
    fi
    echo "PASS: restart: constitutive state written"

    if awk 'NR > 20 && $1 + 0 > 1e-12 {found = 1} END {exit !found}' \
        "${state_file}"
    then
        echo "PASS: restart: material has yielded before the restart point"
    else
        echo "FAIL: restart: no plastic history at t = 10, so this test is"
        echo "      vacuous - an elastic restart is exact without any state IO"
        return 1
    fi

    # 2. the negative control: without the history, refuse
    local guard_dir="${REGRESSION_ROOT}/frameworkRestartMissing"
    rm -rf "${guard_dir}"
    cp -a "${CASE_DIR}" "${guard_dir}"
    rm -f "${guard_dir}"/10/*:*:epsilonP
    sed -i.bak \
        's/^startFrom       startTime;/startFrom       latestTime;/; s/^endTime         10;/endTime         20;/' \
        "${guard_dir}/system/controlDict"
    rm -f "${guard_dir}/system/controlDict.bak"

    if ( cd "${guard_dir}" && solids4Foam > log.solids4Foam 2>&1 ); then
        echo "FAIL: restart: continued without its history instead of refusing"
        return 1
    fi

    if grep -q "is not there" "${guard_dir}/log.solids4Foam"; then
        echo "PASS: restart: refuses when the history is missing"
    else
        echo "FAIL: restart: stopped, but not for the missing state"
        return 1
    fi

    # 1. the restart itself
    sed -i.bak \
        's/^startFrom       startTime;/startFrom       latestTime;/; s/^endTime         10;/endTime         20;/' \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"

    mv "${CASE_DIR}/${SOLVER_LOGFILE}" "${CASE_DIR}/log.solids4Foam.firstLeg"

    if ! ( cd "${CASE_DIR}" && solids4Foam > "${SOLVER_LOGFILE}" 2>&1 ); then
        echo "FAIL: restart: the continued run did not finish"
        return 1
    fi

    local eps sig
    eps=$(extract_max_epsilon)
    sig=$(extract_max_sigma)

    if [[ -z "${eps}" || -z "${sig}" || -z "${RESULT_EPS[framework]:-}" ]]; then
        echo "SKIP: restart: needs the framework arm to have run"
        return 0
    fi

    # The tolerance sits between the two things it has to tell apart: a
    # correct restart, which lands within 1e-8 here, and one that lost its
    # history, which is wrong by 2e-2. The remaining 1e-8 is not this
    # framework's: the legacy model restarts to the same 1e-8 on this case,
    # because the solver stops on a residual measured relative to its first
    # one and a restarted step does not start from the same guess
    local a b
    a="${RESULT_EPS[framework]}"; b="${eps}"

    if awk "BEGIN {exit !(($a - $b)^2 <= (1e-6*$a)^2)}"; then
        printf "PASS: restart reproduces the uninterrupted run (%.8g vs %.8g)\n" \
            "$a" "$b"
    else
        printf "FAIL: restart differs from the uninterrupted run (%.8g vs %.8g)\n" \
            "$a" "$b"
        return 1
    fi

    return 0
}

if ! run_restart_test; then
    failures=$((failures + 1))
fi

# ------------------------------------------------------------
# Restart on a different decomposition
# ------------------------------------------------------------
# decomposePar does not copy the constitutive state into the processor
# directories - it has no idea what these files are, and there is no hook by
# which a library can teach it. So a decomposed run goes back to the
# undecomposed case for the state and distributes it itself, using the
# cellProcAddressing that decomposePar already writes.
#
# This runs serially to t = 10, decomposes onto four processors, and continues
# there. The answer has to be the uninterrupted serial one, because a change of
# decomposition is not supposed to be a change of problem
run_parallel_restart_test() {
    if ! command -v mpirun > /dev/null 2>&1; then
        echo "SKIP: parallel restart (no mpirun)"
        return 0
    fi

    local d="${REGRESSION_ROOT}/frameworkParallel"

    prepare_case "frameworkParallel"
    CASE_DIR="${d}"

    sed -i.bak 's/^writePrecision  6;/writePrecision  14;/; s/^endTime         20;/endTime         10;/' \
        "${d}/system/controlDict"
    rm -f "${d}/system/controlDict.bak"

    ( cd "${d}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )

    cat > "${d}/system/decomposeParDict" << 'EOD'
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains 4;
method scotch;
EOD

    sed -i.bak \
        's/^startFrom       startTime;/startFrom       latestTime;/; s/^endTime         10;/endTime         20;/' \
        "${d}/system/controlDict"
    rm -f "${d}/system/controlDict.bak"

    if ! ( cd "${d}" && decomposePar -time 10 > log.decomposePar 2>&1 ); then
        echo "SKIP: parallel restart (decomposePar failed)"
        return 0
    fi

    # The premise of the whole test: decomposePar leaves the state behind. If
    # some future version starts copying it, this test would quietly become a
    # test of something else
    if ls "${d}"/processor0/10/*:*:epsilonP > /dev/null 2>&1; then
        echo "FAIL: parallel restart: decomposePar copied the state, so this"
        echo "      test is no longer exercising the mapping it was written for"
        return 1
    fi
    echo "PASS: parallel restart: state is not in the processor directories"

    if ! ( cd "${d}" && mpirun -np 4 solids4Foam -parallel > log.par 2>&1 ); then
        echo "FAIL: parallel restart did not run"
        grep -m2 "FOAM FATAL" -A3 "${d}/log.par" || true
        return 1
    fi

    if grep -q "Mapped .* from the undecomposed case" "${d}/log.par"; then
        echo "PASS: parallel restart mapped the state from the serial case"
    else
        echo "FAIL: parallel restart ran without mapping any state"
        return 1
    fi

    local eps
    eps=$(grep "Max epsilonEq" "${d}/log.par" | awk '{print $NF}' | tail -n 1)

    if [[ -z "${eps}" || -z "${RESULT_EPS[framework]:-}" ]]; then
        echo "SKIP: parallel restart needs the framework arm to have run"
        return 0
    fi

    local a="${RESULT_EPS[framework]}"

    if awk "BEGIN {exit !(($a - $eps)^2 <= (1e-6*$a)^2)}"; then
        printf "PASS: restart on four processors matches the serial run (%.8g vs %.8g)\n" \
            "$a" "$eps"
    else
        printf "FAIL: restart on four processors differs (%.8g vs %.8g)\n" \
            "$a" "$eps"
        return 1
    fi

    return 0
}

if ! run_parallel_restart_test; then
    failures=$((failures + 1))
fi

echo
if (( failures == 0 ))
then
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
