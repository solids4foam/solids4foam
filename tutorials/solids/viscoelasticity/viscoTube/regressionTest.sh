#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

# For solids4Foam::foamFlavour, which picks the legacy reference for this fork
source "${SOLIDS4FOAM_SCRIPTS}"

# viscousHookeanElastic is history dependent through its Maxwell arms, and
# this is the only end-to-end comparison of the mechanicalConstitutiveLaw
# framework's law against the removed legacy one. There is no pressure
# smoothing here, so the framework must reproduce the legacy answer
APPROACHES=(
    main
)

# ============================================================
# viscoTube regression test
# Checks order of magnitude of strain and von Mises stress
# ============================================================

# Order-of-magnitude tolerances
EPS_MIN=1e-5
EPS_MAX=1e-3

SIGMA_MIN=1e6
SIGMA_MAX=1e8

# The final extrema of the removed legacy mechanicalModel, logged to fourteen
# figures, from the last commit that had it (mcl-stage8-coverage, c3a92b3d),
# per fork. The framework matched them to the eight digits compared, and is
# held to the 1e-6 relative that comparison used. The Maxwell arms are history,
# so a relaxation error accumulates over the run and shows here in a way the
# unit checks, which compare two trial states from one rest state, cannot see
case "$(solids4Foam::foamFlavour)" in
    com)
        LEGACY_EPS=0.00019987366
        LEGACY_SIG=11480085
        ;;
    org)
        LEGACY_EPS=0.00019987366
        LEGACY_SIG=11480085
        ;;
    foamextend)
        LEGACY_EPS=0.00020025214
        LEGACY_SIG=11501844
        ;;
esac
LEGACY_REL_TOL=1e-6

# Log files
SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
CONSTITUTIVE_LOGFILE="log.Test-mechanicalConstitutiveLaw"

echo "============================================================"
echo "viscoTube regression test"
echo "epsilonEq expected: ${EPS_MIN} < eps < ${EPS_MAX}"
echo "sigmaEq expected  : ${SIGMA_MIN} < sigma < ${SIGMA_MAX}"
echo "============================================================"
echo

# Exercise the framework's own checks on this case. Its law is
# viscousHookeanElastic, the first law whose response depends on the time
# increment, so this is the runtime coverage of the inputs object carrying dt
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
            > "${CONSTITUTIVE_LOGFILE}" 2>&1 )
    then
        local n_passed
        n_passed=$(grep -c 'PASS:' "${CASE_DIR}/${CONSTITUTIVE_LOGFILE}" || true)

        if (( n_passed == 0 )); then
            echo "FAIL: mechanicalConstitutiveLaw checks reported no checks"
            return 1
        fi

        echo "PASS: mechanicalConstitutiveLaw checks (${n_passed} checks)"
        return 0
    fi

    echo "FAIL: mechanicalConstitutiveLaw checks"
    grep 'FAIL:' "${CASE_DIR}/${CONSTITUTIVE_LOGFILE}" || true
    return 1
}

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

    # Every arm at the same precision. Comparing a value logged at six
    # significant figures against one logged at fourteen measures the log
    # format and calls the difference a regression
    sed -i.bak 's/^writePrecision  6;/writePrecision  14;/' \
        "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"

    # restart yes makes the solid model write the kinematic history the
    # restart test below needs; the main arm asks too, so that the two are set
    # up the same
    sed -i.bak \
        's/^    nCorrectors     1000;/    restart yes;\n    nCorrectors     1000;/' \
        "${CASE_DIR}/constant/solidProperties"
    rm -f "${CASE_DIR}/constant/solidProperties.bak"

    if ! grep -q 'restart yes' "${CASE_DIR}/constant/solidProperties"; then
        echo "FAIL: could not set restart in solidProperties"
        exit 1
    fi
}

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_epsilon() {
    grep "Max epsilonEq" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk -F '=' '{print $2}' \
        | tr -d '[:space:]'
}

extract_max_sigma() {
    grep "Max sigmaEq (von Mises stress)" "${CASE_DIR}/${SOLVER_LOGFILE}" \
        | tail -n 1 \
        | awk -F '=' '{print $2}' \
        | tr -d '[:space:]'
}

# ------------------------------------------------------------
# Run and check each approach
# ------------------------------------------------------------

failures=0
declare -A RESULT_E
declare -A RESULT_S

for approach in "${APPROACHES[@]}"; do
    echo
    echo "------------------------------------------------------------"
    echo "Testing approach: ${approach}"
    echo "------------------------------------------------------------"

    prepare_case "${approach}"
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )

    if grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null
    then
        echo "PASS: ${approach} took its material from the framework"
    else
        echo "FAIL: ${approach} constructed no mechanical constitutive law"
        failures=$((failures + 1))
    fi

    epsilon=$(extract_max_epsilon)
    sigma=$(extract_max_sigma)

    if [[ -z "${epsilon}" || -z "${sigma}" ]]; then
        echo "FAIL: ${approach} could not extract epsilonEq or sigmaEq"
        failures=$((failures + 1))
        continue
    fi

    RESULT_E["${approach}"]="${epsilon}"
    RESULT_S["${approach}"]="${sigma}"

    if awk "BEGIN {exit !(${epsilon} > ${EPS_MIN} && ${epsilon} < ${EPS_MAX})}"; then
        printf "PASS: %s Max epsilonEq = %.6g\n" "${approach}" "${epsilon}"
    else
        printf "FAIL: %s Max epsilonEq = %.6g\n" "${approach}" "${epsilon}"
        failures=$((failures + 1))
    fi

    if awk "BEGIN {exit !(${sigma} > ${SIGMA_MIN} && ${sigma} < ${SIGMA_MAX})}"; then
        printf "PASS: %s Max sigmaEq = %.6g\n" "${approach}" "${sigma}"
    else
        printf "FAIL: %s Max sigmaEq = %.6g\n" "${approach}" "${sigma}"
        failures=$((failures + 1))
    fi

    if ! run_constitutive_test; then
        failures=$((failures + 1))
    fi

    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
done

if [[ -n "${RESULT_E[main]:-}" ]]; then
    for q in eps sig; do
        if [[ "${q}" == "eps" ]]; then
            a="${LEGACY_EPS}"; b="${RESULT_E[main]}"; n="epsilonEq"
        else
            a="${LEGACY_SIG}"; b="${RESULT_S[main]}"; n="sigmaEq"
        fi

        if awk "BEGIN {exit !(($a - $b)^2 <= (${LEGACY_REL_TOL}*$a)^2)}"; then
            printf "PASS: %s matches the legacy model (%.8g vs %.8g)\n" "$n" "$b" "$a"
        else
            printf "FAIL: %s differs from the legacy model (%.8g vs %.8g)\n" "$n" "$b" "$a"
            failures=$((failures + 1))
        fi
    done
else
    echo "FAIL: the main arm produced nothing to compare with the legacy model"
    failures=$((failures + 1))
fi


# ------------------------------------------------------------
# Restart
# ------------------------------------------------------------
# This law's history is one internal stress per Maxwell arm, and how many arms
# there are comes from the material's own dictionary - so it is the case that
# checks a law can declare a variable number of state fields and get all of
# them back. A viscous material is also the clearest possible test of whether
# history survived: relaxation is the whole behaviour, and a run that forgot
# how far each arm had relaxed starts again from an unstressed state
run_restart_test() {
    local d="${REGRESSION_ROOT}/frameworkRestart"

    prepare_case "frameworkRestart"
    CASE_DIR="${d}"

    sed -i.bak 's/^endTime         7000;/endTime         3500;/' \
        "${d}/system/controlDict"
    rm -f "${d}/system/controlDict.bak"

    ( cd "${d}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
        echo "FAIL: restart: the first leg did not run"
        return 1
    }

    # One file per arm, plus the relaxing deviatoric stress
    local nArms
    nArms=$(ls "${d}"/3500/*IntegrationPointTopology_h[0-9]* 2>/dev/null | wc -l)

    if ! ls "${d}"/3500/*IntegrationPointTopology_s > /dev/null 2>&1 || (( nArms == 0 )); then
        echo "FAIL: restart: the viscous history was not written"
        return 1
    fi
    echo "PASS: restart: the viscous history is written (${nArms} arm(s))"

    sed -i.bak \
        's/^startFrom       startTime;/startFrom       latestTime;/; s/^endTime         3500;/endTime         7000;/' \
        "${d}/system/controlDict"
    rm -f "${d}/system/controlDict.bak"

    mv "${d}/${SOLVER_LOGFILE}" "${d}/log.solids4Foam.firstLeg"

    if ! ( cd "${d}" && solids4Foam > "${SOLVER_LOGFILE}" 2>&1 ); then
        echo "FAIL: restart: the continued run did not finish"
        grep -m1 "FOAM FATAL" -A4 "${d}/${SOLVER_LOGFILE}" || true
        return 1
    fi

    local eps sig
    eps=$(extract_max_epsilon)
    sig=$(extract_max_sigma)

    if [[ -z "${eps}" || -z "${RESULT_S[main]:-}" ]]; then
        echo "SKIP: restart needs the main arm to have run"
        return 0
    fi

    # Stress, not strain: strain is driven by the load and comes back whatever
    # the history did, while the stress is what the relaxation determines
    local a="${RESULT_S[main]}"

    if awk "BEGIN {exit !(($a - $sig)^2 <= (1e-6*$a)^2)}"; then
        printf "PASS: restart reproduces the uninterrupted run (%.8g vs %.8g)\n" \
            "$a" "$sig"
    else
        printf "FAIL: restart differs from the uninterrupted run (%.8g vs %.8g)\n" \
            "$a" "$sig"
        return 1
    fi

    return 0
}

if ! run_restart_test; then
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
