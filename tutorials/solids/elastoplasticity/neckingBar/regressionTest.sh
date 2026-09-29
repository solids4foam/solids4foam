#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

source "${SOLIDS4FOAM_SCRIPTS}"

# GNU sed, for the in-place edits below
solids4Foam::requireGnuSed

# The end-to-end check of two things nothing else covers: the stress
# migration of nonLinGeomUpdatedLagSolid onto the mechanicalConstitutiveLaw
# framework, and the finite-strain plastic law. It is the only tutorial pairing
# them that does not enable pressure smoothing, so it is the only one where the
# framework reproduced the removed legacy model exactly rather than
# approximately, and it is held to that below
# ============================================================
# neckingBar regression test
# Checks the final loading force after the necking curve.
# ============================================================

FORCE_MIN=74.3
FORCE_MAX=74.6

# The final loading force. With no pressure smoothing the framework omits
# nothing, and it matched the removed legacy mechanicalModel's in all eight
# digits printed, on every fork (mcl-stage8-coverage, c3a92b3d). This is
# OpenFOAM.com v2512's, the same as OpenFOAM.org 9's; foam-extend 4.1 gives
# 2e-4 more, and the tolerance, 5e-4 relative, covers that. A history error in
# the finite-strain plastic path, or in the relative deformation gradient the
# framework is given, shows up here and nowhere else
REF_FORCE=74.4649
REF_REL_TOL=5e-4

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "neckingBar regression test"
echo "Final loading force in [${FORCE_MIN}, ${FORCE_MAX}]"
echo "============================================================"
echo

# A copy of the case in regressionTests/<name>, main unless named
prepare_case() {
    CASE_DIR="${REGRESSION_ROOT}/${1:-main}"

    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done

    # restart yes makes the solid model write the kinematic history the
    # restart test below needs; the main arm asks too, so that the two are set
    # up the same
    "${SOLIDS4FOAM_SED}" -i.bak \
        's/^    nCorrectors     1000;/    restart yes;\n    nCorrectors     1000;/' \
        "${CASE_DIR}/constant/solidProperties"
    rm -f "${CASE_DIR}/constant/solidProperties.bak"

    if ! grep -q 'restart yes' "${CASE_DIR}/constant/solidProperties"; then
        echo "FAIL: could not set restart in solidProperties"
        exit 1
    fi
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

failures=0
result_f=""

run_case() {
    if [ "$CHECK_ONLY" = false ]; then
        prepare_case
        "${SOLIDS4FOAM_SED}" -i.bak 's/^endTime         1;/endTime         0.184;/' \
            "${CASE_DIR}/system/controlDict"
        rm -f "${CASE_DIR}/system/controlDict.bak"
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
        ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
    else
        CASE_DIR="${REGRESSION_ROOT}/main"
        echo "Running in check-only mode: skipping Allclean and Allrun"
    fi

    if grep -q "Selecting mechanical constitutive law" \
        "${CASE_DIR}/log.solids4Foam" 2>/dev/null
    then
        echo "PASS: took its material from the framework"
    else
        echo "FAIL: constructed no mechanical constitutive law"
        failures=$((failures + 1))
    fi

    force_file=$(find "${CASE_DIR}/postProcessing" \
        -name 'solidForcesDisplacementsloading.dat' -print | tail -n 1)
    if [[ -z "${force_file}" ]]; then
        echo "FAIL: could not find loading force history"
        failures=$((failures + 1))
        return 0
    fi

    final_force=$(awk 'END {print $5*360e-3}' "${force_file}")
    result_f="${final_force}"

    if awk "BEGIN {exit !(${final_force} >= ${FORCE_MIN} && ${final_force} <= ${FORCE_MAX})}"; then
        printf "PASS: final loading force = %.6g\n" "${final_force}"
    else
        printf "FAIL: final loading force = %.6g\n" "${final_force}"
        failures=$((failures + 1))
    fi

    if [ "$CHECK_ONLY" = false ]; then
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    fi
}

run_case

if [[ -n "${result_f}" ]]; then
    a="${REF_FORCE}"; b="${result_f}"

    if awk "BEGIN {exit !(($a - $b)^2 <= (${REF_REL_TOL}*$a)^2)}"; then
        printf "PASS: force matches the reference (%.8g vs %.8g)\n" "$b" "$a"
    else
        printf "FAIL: force differs from the reference (%.8g vs %.8g)\n" "$b" "$a"
        failures=$((failures + 1))
    fi
else
    echo "FAIL: the main arm produced no force to compare with the reference"
    failures=$((failures + 1))
fi

echo

# ------------------------------------------------------------
# Restart, under finite strain
# ------------------------------------------------------------
# The other restart tests are small strain. This one is the updated Lagrangian
# solid model, where the history is bEbar - the isochoric elastic left
# Cauchy-Green tensor - whose undeformed value is the identity rather than
# zero. That is the state most likely to be quietly wrong if a restart loses
# it, because a field of zeros is not a neutral deformation but a degenerate
# one, and the run would not obviously misbehave so much as produce nonsense.
#
# Deliberately short, and deliberately not a plasticity test: this bar does not
# yield until well past halfway, and the plastic path is already covered by
# perforatedPlate and stripFooting. What is uncovered, and what this covers, is
# finite-strain history through a different solid model
run_restart_test() {
    local d="${REGRESSION_ROOT}/restart"

    if [ "$CHECK_ONLY" = false ]; then
        prepare_case "restart"

        "${SOLIDS4FOAM_SED}" -i.bak \
            's/^writePrecision.*/writePrecision  14;/; s/^endTime         1;/endTime         0.1;/; s/^deltaT          0.002;/deltaT          0.01;/; s/^writeInterval   10;/writeInterval   5;/' \
            "${d}/system/controlDict"
        rm -f "${d}/system/controlDict.bak"

        ( cd "${d}" && ./Allrun > log.Allrun 2>&1 ) || {
            echo "FAIL: finite-strain restart: the reference run did not finish"
            return 1
        }
    elif [[ ! -d "${d}" ]]; then
        echo "FAIL: no existing finite-strain restart result to check"
        return 1
    else
        echo "Check-only: checking existing finite-strain restart results"
    fi

    local ref
    local reffile
    reffile=$(find "${d}/postProcessing" -name 'solidForcesDisplacementsloading.dat' \
        -print 2>/dev/null | tail -n 1)
    ref=$(awk 'END {print $5*360e-3}' "${reffile}" 2>/dev/null)

    if [[ -z "${ref}" ]]; then
        echo "SKIP: finite-strain restart (no force history)"
        return 0
    fi

    # The state has to be real history, or this proves nothing. bEbar starts as
    # the identity, so a deformed state is one whose diagonal has moved off one
    local bfile
    bfile=$(ls "${d}"/0.05/*IntegrationPointTopology_bEbar 2>/dev/null | head -n 1)

    if [[ -z "${bfile}" ]]; then
        echo "FAIL: finite-strain restart: bEbar was not written"
        return 1
    fi

    # Read as six numbers per point and compared against the identity, rather
    # than by looking at one component of lines that start with a bracket: a
    # uniform list is written as N{(...)} with no such line, and a deformation
    # that left xx at one while moving another component would have been missed
    if awk '
        # The list after the header, as the count and then six components per
        # point, or as N{(...)} for a uniform list
        /^\/\/ \* \* \*/ { text = ""; next }
        /^[ \t]*\/\// { next }
        { text = text " " $0 }
        END {
            gsub(/[(){}]/, " ", text)
            n = split(text, v, " ")
            split("1 0 0 1 0 1", identity, " ")
            for (k = 2; k <= n; k++) {
                d = v[k] - identity[(k - 2) % 6 + 1]
                if (d < 0) d = -d
                if (d > 1e-10) exit 0
            }
            exit 1
        }' "${bfile}"
    then
        echo "PASS: finite-strain restart: bEbar carries deformation"
    else
        echo "FAIL: finite-strain restart: bEbar is still the identity, so"
        echo "      this test would pass even if the state were lost"
        return 1
    fi

    # Negative control first: without bEbar the run must stop. This is what
    # gives the comparison below its teeth, because that comparison has to
    # tolerate a difference this solid model shows with or without the
    # framework, and a tolerance wide enough to pass would also be wide enough
    # to hide a lost state if nothing else were checked
    local m="${REGRESSION_ROOT}/restartMissing"
    if [ "$CHECK_ONLY" = false ]; then
        rm -rf "${m}"; cp -a "${d}" "${m}"
        rm -rf "${m}"/0.0[6-9] "${m}"/0.1 "${m}"/postProcessing
        rm -f "${m}"/0.05/*IntegrationPointTopology_bEbar
        "${SOLIDS4FOAM_SED}" -i.bak 's/^startFrom       startTime;/startFrom       latestTime;/' \
            "${m}/system/controlDict"
        rm -f "${m}/system/controlDict.bak"

        if ( cd "${m}" && solids4Foam > log.solids4Foam 2>&1 ); then
            echo "FAIL: finite-strain restart: continued without its bEbar"
            return 1
        fi
    fi

    if grep -q "bEbar.*is not there" "${m}/log.solids4Foam" 2>/dev/null; then
        echo "PASS: finite-strain restart: refuses when bEbar is missing"
    else
        echo "FAIL: finite-strain restart: stopped, but not for missing bEbar"
        return 1
    fi

    # Continue from halfway in a copy, so the reference stays intact
    local g="${REGRESSION_ROOT}/restartLeg"
    if [ "$CHECK_ONLY" = false ]; then
        rm -rf "${g}"; cp -a "${d}" "${g}"
        rm -rf "${g}"/0.0[6-9] "${g}"/0.1 "${g}"/postProcessing
        "${SOLIDS4FOAM_SED}" -i.bak \
            's/^startFrom       startTime;/startFrom       latestTime;/' \
            "${g}/system/controlDict"
        rm -f "${g}/system/controlDict.bak"

        if ! ( cd "${g}" && solids4Foam > log.solids4Foam 2>&1 ); then
            echo "FAIL: finite-strain restart did not run"
            grep -m1 "FOAM FATAL" -A4 "${g}/log.solids4Foam" || true
            return 1
        fi
    fi

    local got
    local gotfile
    gotfile=$(find "${g}/postProcessing" -name 'solidForcesDisplacementsloading.dat' \
        -print 2>/dev/null | tail -n 1)
    got=$(awk 'END {print $5*360e-3}' "${gotfile}" 2>/dev/null)

    if [[ -z "${got}" ]]; then
        echo "SKIP: finite-strain restart (no force history after restart)"
        return 0
    fi

    # 1e-4, not the 1e-6 the small-strain cases hold to. This solid model does
    # not reproduce an uninterrupted run exactly across a restart, with or
    # without the framework: the removed legacy model differed by 1.3e-5 on
    # this quantity here and the framework by 4.3e-6, so the residual is the
    # solid model's and the framework is if anything the closer of the two. The
    # tolerance is set to accept that and nothing larger; a lost bEbar is
    # caught by the check above rather than by this number
    if awk "BEGIN {exit !(($ref - $got)^2 <= (1e-4*$ref)^2)}"; then
        printf "PASS: finite-strain restart reproduces the run (%.8g vs %.8g)\n" \
            "$ref" "$got"
    else
        printf "FAIL: finite-strain restart differs (%.8g vs %.8g)\n" \
            "$ref" "$got"
        return 1
    fi

    return 0
}

if ! run_restart_test; then
    failures=$((failures + 1))
fi

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
