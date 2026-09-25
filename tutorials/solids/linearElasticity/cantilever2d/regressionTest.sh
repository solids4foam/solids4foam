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

# ------------------------------------------------------------
# The answers of the removed legacy mechanicalModel
# ------------------------------------------------------------
# From the last commit that had it (mcl-stage8-coverage, c3a92b3d), per fork.
# Each is the legacy side of a comparison the framework passed there, and is
# held to the tolerance that comparison used. Fields are recorded as the max
# and mean component magnitude of their internal values as written
case "$(solids4Foam::foamFlavour)" in
    com|org)
        LEGACY_POINT_D_MAX=0.015999512315885
        LEGACY_POINT_D_MEAN=0.00207382900687327
        LEGACY_EXPLICIT_POINT_D_MAX=4.50000000005551e-06
        LEGACY_EXPLICIT_POINT_D_MEAN=2.69525341068081e-09
        LEGACY_EXPLICIT_DELTAT_LINE="Setting deltaT = 8.830682704707434e-08, maxCo = 0.1"
        # unsCoupled runs on foam-extend only
        LEGACY_UNSCOUPLED_EPS=""
        ;;
    foamextend)
        LEGACY_POINT_D_MAX=0.0159995123112567
        LEGACY_POINT_D_MEAN=0.00207382900625648
        LEGACY_EXPLICIT_POINT_D_MAX=4.5e-06
        LEGACY_EXPLICIT_POINT_D_MEAN=2.6952526893628e-09
        LEGACY_EXPLICIT_DELTAT_LINE="Setting deltaT = 8.830682704706176e-08, maxCo = 0.1"
        LEGACY_UNSCOUPLED_EPS=0.00044921
        ;;
esac

# The same on every fork
declare -A LEGACY_SIGMA0_D_DIGEST=(
    [dict]=1c2bf5df020259760a2afc4badb6cf86e101c1e5
    [field]=2c8072a41810075abbdbd5912646a309423536bf
    [both]=1c2bf5df020259760a2afc4badb6cf86e101c1e5
)

# unsCoupled: the final max epsilonEq. The same problem with the same material
# constants read two ways, so the two agreed far more closely than the band
# allows, in every figure logged. The log gives five, so a round-off difference
# on another machine can move the last one; the tolerance, 1e-4 relative, is
# a few units in that figure
UNSCOUPLED_REL_TOL=1e-4

# The vertex-centred model's final pointD, implicit and explicit: the same
# linear elastic problem, so the framework reproduced the legacy field to
# round-off, and exactly as written. CI measures 3.6e-10 relative against the
# values recorded on macOS, on every fork. The tolerance, 1e-6 of the largest
# component, allows for that, and is ten times below the 1e-5 that a 0.001%
# change in E makes
VERTEX_CENTRED_DISP_REL_TOL=1e-6

# sigma0, given in the law's dictionary, as a field, and as both: the final D
# agreed with the legacy model's exactly, in all six figures written, and is
# recorded as a digest of its internal values

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

APPROACHES=(
    petscSnes
    unsCoupled
    highOrder-movingLeastSquares
    highOrder-kExactLeastSquares
    highOrderJacobian
    vertexCentred
)

echo "============================================================"
echo "cantilever2d regression test"
echo "Max epsilonEq in [${EPS_MIN}, ${EPS_MAX}]"
echo "Max sigmaEq   in [${SIGMA_MIN}, ${SIGMA_MAX}]"
echo "High-order DDifference LInf < ${HIGH_ORDER_DISP_TOL}"
echo "Against the legacy model: unsCoupled max epsilonEq to ${UNSCOUPLED_REL_TOL},"
echo "  vertex-centred pointD to ${VERTEX_CENTRED_DISP_REL_TOL}, sigma0 D exactly"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    # A result kept from the vertex-centred arm of an earlier run would let
    # this run's comparison pass without the arm having produced anything
    rm -f "${REGRESSION_ROOT}"/pointD.vertexCentred*

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

# Run this case with an initial stress, and require the answer of the removed
# legacy mechanical law exactly.
#
# sigma0 is the first prescribed state the framework declares. It reaches the
# law by two routes - a uniform value in the law's dictionary, and a field the
# case supplies - and both are checked here because they are read by different
# code. A non-uniform field is used for the second so that a wrong
# cell-to-integration-point map changes the answer, which a uniform field
# would hide.
#
# This case is the host because it converges tightly enough for the framework
# and the legacy model to have agreed to the last digit, which makes the check
# a statement about the model rather than about a solver tolerance
run_sigma0_check() {
    local mode="$1"
    local d="${REGRESSION_ROOT}/sigma0-${mode}"

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

        if ! grep -q "sigma0 sigma0" "${d}/constant/mechanicalProperties"
        then
            echo "FAIL: could not give sigma0 in the law's dictionary"
            return 1
        fi
    fi

    (
        cd "${d}" || exit 1

        solids4Foam::convertCaseFormat . > log.convert 2>&1

        blockMesh > log.blockMesh 2>&1 || exit 1

        if [[ "${mode}" == "field" || "${mode}" == "both" ]]; then
            write_sigma0_field . || exit 1
        fi

        solids4Foam > log.solids4Foam 2>&1 || exit 1
    ) || { echo "FAIL: sigma0 ${mode} could not run"; return 1; }

    if [[ "${mode}" == "dict" || "${mode}" == "both" ]]; then
        if ! grep -q "Uniform initial stress sigma0" "${d}/log.solids4Foam"
        then
            echo "FAIL: sigma0 ${mode}: the law did not read sigma0 from the dict"
            return 1
        fi
    else
        if ! grep -q "Prescribed state 'sigma0'" "${d}/log.solids4Foam"
        then
            echo "FAIL: sigma0 ${mode}: the law did not read the sigma0 field"
            return 1
        fi
    fi

    local t
    t=$(solids4Foam::latestTime "${d}")

    if [[ -z "${t}" || ! -f "${d}/${t}/D" ]]; then
        echo "FAIL: sigma0 ${mode} produced no D field"
        return 1
    fi

    # sigma0 must actually have changed the answer, or agreement is vacuous
    if diff -q "${d}/${t}/D" "${SIGMA0_BASELINE_D}" > /dev/null 2>&1
    then
        echo "FAIL: sigma0 ${mode} left the solution unchanged"
        return 1
    fi

    # The legacy model read a sigma0 field and then assigned any dictionary
    # sigma0 over the whole of it, so when a case carried both, the dictionary
    # was what took effect. The framework must resolve the two the same way
    # round
    if [[ "${mode}" == "dict" ]]; then
        SIGMA0_DICT_D="${d}/${t}/D"
    elif [[ "${mode}" == "both" ]]; then
        if ! diff -q "${d}/${t}/D" "${SIGMA0_DICT_D}" > /dev/null
        then
            echo "FAIL: sigma0 both: the field was not overridden by the dict"
            return 1
        fi
    fi

    local digest
    digest=$(internal_field_digest "${d}/${t}/D" || true)

    if [[ -z "${digest}" \
        || "${digest}" != "${LEGACY_SIGMA0_D_DIGEST[${mode}]}" ]]
    then
        echo "FAIL: sigma0 ${mode}: D differs from the legacy model's" \
            "(digest ${digest:-none}, legacy ${LEGACY_SIGMA0_D_DIGEST[${mode}]})"
        return 1
    fi

    echo "PASS: sigma0 ${mode}: D matches the legacy model's exactly"
    return 0
}

# A prescribed field has to survive a restart. It is written into 0 and the
# case then restarts from a later time whose directory does not contain it, so
# this is the check that the field is looked for where it actually lives.
#
# It compares the framework against itself, restarted against continuous,
# rather than against the legacy model. That looked for sigma0 only beside the
# fields it restarted from and so lost it, keeping it only when an earlier run
# had already written it forward. That was a hole rather than a behaviour worth
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

        # restart no, deliberately. These arms restart, and a case that
        # restarts without saying anything about it is refused - but this law
        # is written in total strain and genuinely does not need the kinematic
        # history, so saying no is the honest answer rather than a way round
        # the check. 'restart yes' would also make legacy's sigma0 be written
        # into the time this resumes from, which is the very thing the check
        # below requires to be absent
        sed -i \
            's|^\( *\)nCorrectors|\1restart no;\n\1nCorrectors|' \
            "${root}/${d}/constant/solidProperties"

        if ! grep -q "restart no;" "${root}/${d}/constant/solidProperties"
        then
            echo "FAIL: sigma0 restart check could not set restart in ${d}"
            return 1
        fi
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

    for d in restartContinuous restartRestarted restartNone; do
        if ! grep -q "Selecting mechanical constitutive law" \
            "${root}/${d}/log.solids4Foam"
        then
            echo "FAIL: sigma0 restart arm ${d} constructed no mechanical" \
                "constitutive law"
            return 1
        fi
    done

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
    t=$(solids4Foam::latestTime "${d}")
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

# The case directory is cleaned between approaches, so the unsCoupled result
# is captured as it goes and compared with the legacy answer afterwards
UNSCOUPLED_EPS=""

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

    if [[ "${approach}" == unsCoupled ]]; then
        UNSCOUPLED_EPS="${epsilon}"
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

# A digest of a field's internal values as written, independent of the file
# header, which names the OpenFOAM version that wrote it
internal_field_digest() {
    python3 - "$1" << 'PYEOF'
import hashlib
import re
import sys

text = open(sys.argv[1]).read()
field = re.search(r"\binternalField\s+(.*?);\s*\n\s*boundaryField", text, re.DOTALL)
if not field:
    sys.exit(f"cannot find the internalField in {sys.argv[1]}")
print(hashlib.sha1(" ".join(field.group(1).split()).encode()).hexdigest())
PYEOF
}

# The largest magnitude of any component of a field's internal values, and the
# mean magnitude, as "max<TAB>mean"
internal_field_norms() {
    python3 - "$1" << 'PYEOF'
import re
import sys

number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
text = open(sys.argv[1]).read()

uniform = re.search(
    r"\binternalField\s+uniform\s+(\([^)]*\)|" + number + r")\s*;", text
)
if uniform:
    body = uniform.group(1)
else:
    field = re.search(
        r"\binternalField\s+nonuniform\s+List<\w+>\s+\d+\s*\((.*?)\n\)\s*;",
        text,
        re.DOTALL,
    )
    if not field:
        sys.exit(f"cannot parse internalField in {sys.argv[1]}")
    body = field.group(1)

values = [abs(float(x)) for x in re.findall(number, body)]
if not values:
    sys.exit(f"empty internalField in {sys.argv[1]}")
print(f"{max(values):.15g}\t{sum(values)/len(values):.15g}")
PYEOF
}

# A field against the removed legacy model's, through the norms above. Both
# differences are bounded by the largest pointwise difference, so a field that
# agrees with the legacy one to tol times its largest value passes, and one
# that does not is caught by at least one of the two in all but contrived cases
check_field_against_legacy() {
    local label="$1"
    local file="$2"
    local legacy_max="$3"
    local legacy_mean="$4"
    local tol="$5"
    local norms field_max field_mean

    if [[ ! -f "${file}" ]] || ! norms=$(internal_field_norms "${file}"); then
        echo "FAIL: ${label}: no field to compare with the legacy model"
        return 1
    fi

    read -r field_max field_mean <<< "${norms}"

    if awk "BEGIN {
            a = ${field_max} - ${legacy_max}; if (a < 0) a = -a
            b = ${field_mean} - ${legacy_mean}; if (b < 0) b = -b
            exit !(${field_max} > 0 && a <= ${tol}*${legacy_max} \
                && b <= ${tol}*${legacy_max})
        }"
    then
        printf "PASS: %s matches the legacy model: max %.15g (%.15g), mean %.15g (%.15g)\n" \
            "${label}" "${field_max}" "${legacy_max}" "${field_mean}" "${legacy_mean}"
        return 0
    fi

    printf "FAIL: %s differs from the legacy model: max %.15g (%.15g), mean %.15g (%.15g), tolerance %s\n" \
        "${label}" "${field_max}" "${legacy_max}" "${field_mean}" "${legacy_mean}" "${tol}"
    return 1
}

# The explicit path of the vertex-centred model, which no tutorial runs:
# twenty steps of the cantilever from rest, in a case directory of its own.
# The time step comes from the wave speed, so it is also a check that the
# framework gives the density and stiffness the legacy model did
run_vertex_centred_explicit_check() {
    local d="${REGRESSION_ROOT}/vertexCentredExplicit"
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

    sed -i \
        "s|^SOLIDS4FOAM_ROOT := .*|SOLIDS4FOAM_ROOT := ${SOLIDS4FOAM_ROOT_ABS}|" \
        "${d}/src/Make/options"

    sed -i \
        's|solutionAlgorithm PETScSNES;|solutionAlgorithm explicit;|' \
        "${d}/constant/solidProperties.vertexCentred"

    # Twenty steps, and only the last one written. The time step is set
    # by the model, so the run is stopped by step count rather than time
    sed -i \
        -e 's|^stopAt .*|stopAt          nextWrite;|' \
        -e 's|^deltaT .*|deltaT          1e-8;|' \
        -e 's|^writeControl .*|writeControl    timeStep;|' \
        -e 's|^writeInterval .*|writeInterval   20;|' \
        -e 's|^writePrecision .*|writePrecision  16;|' \
        "${d}/system/controlDict"

    ( cd "${d}" && ./Allrun vertexCentred > "${ALLRUN_LOGFILE}" 2>&1 ) \
        || { echo "FAIL: the vertexCentred explicit run failed"; return 1; }

    if solids4Foam::regressionCaseSkipped "${d}/${ALLRUN_LOGFILE}"
    then
        echo "Skipping vertexCentred explicit because it is unavailable here"
        return 0
    fi

    if ! grep -q "Selecting mechanical constitutive law" \
        "${d}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: vertexCentred explicit constructed no mechanical" \
            "constitutive law"
        return 1
    fi

    local deltaT_line
    deltaT_line=$(grep "Setting deltaT" "${d}/${SOLVER_LOGFILE}" || true)

    if [[ -z "${deltaT_line}" ]]; then
        echo "FAIL: vertexCentred explicit did not run explicitly"
        return 1
    fi

    if [[ "${deltaT_line}" != "${LEGACY_EXPLICIT_DELTAT_LINE}" ]]; then
        echo "FAIL: vertexCentred explicit chose a different time step from" \
            "the legacy model ('${deltaT_line}', legacy" \
            "'${LEGACY_EXPLICIT_DELTAT_LINE}')"
        return 1
    fi
    echo "PASS: vertexCentred explicit chose the legacy model's time step"

    local t
    t=$(solids4Foam::latestTime "${d}")

    check_field_against_legacy "vertexCentred explicit pointD" \
        "${d}/${t}/pointD" \
        "${LEGACY_EXPLICIT_POINT_D_MAX}" "${LEGACY_EXPLICIT_POINT_D_MEAN}" \
        "${VERTEX_CENTRED_DISP_REL_TOL}"
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

        # The vertex-centred arm is compared with the legacy answer, which the
        # default six significant figures would only do to six figures
        if [[ "${RUN_APPROACH}" == vertexCentred ]]; then
            sed -i 's|^writePrecision .*|writePrecision  16;|' \
                "${CASE_DIR}/system/controlDict"
        fi

        ( cd "${CASE_DIR}" && ./Allrun "${RUN_APPROACH}" > "${ALLRUN_LOGFILE}" 2>&1 )

        if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
            echo "Skipping ${approach} because it is unavailable in this environment"
            continue
        fi

        # Recorded so that a vertex-centred arm which ran but wrote no
        # displacement fails below, rather than quietly switching the
        # comparison off
        if [[ "${RUN_APPROACH}" == vertexCentred ]]; then
            VERTEX_CENTRED_RAN=true
        fi

        if [[ -f "${CASE_DIR}/${SOLVER_LOGFILE}" ]] \
            && ! grep -q "Selecting mechanical constitutive law" \
                "${CASE_DIR}/${SOLVER_LOGFILE}"
        then
            echo "FAIL: ${approach} constructed no mechanical constitutive law"
            failures=$((failures + 1))
        fi

        if ! check_solver_extrema "${approach}"; then
            failures=$((failures + 1))
        fi

        if [[ "${RUN_APPROACH}" == highOrder* ]] && ! check_high_order_errors; then
            failures=$((failures + 1))
        fi

        # Keep the vertex-centred result, since the case directory is
        # cleaned before the next arm runs
        if [[ "${RUN_APPROACH}" == vertexCentred ]]; then
            t=$(solids4Foam::latestTime "${CASE_DIR}")
            if [[ -n "${t}" && -f "${CASE_DIR}/${t}/pointD" ]]; then
                cp "${CASE_DIR}/${t}/pointD" \
                    "${REGRESSION_ROOT}/pointD.${approach}"
            fi
        fi

    done

    # As for unsCoupled, but on the field rather than on one extremum, since
    # the framework reproduced the legacy field to round-off everywhere.
    # Scheduled on whether the arm ran, not on whether it wrote output: an arm
    # that ran and wrote no pointD is a failure, not a reason to skip
    if [[ -n "${VERTEX_CENTRED_RAN:-}" ]]
    then
        if ! check_field_against_legacy "vertexCentred pointD" \
            "${REGRESSION_ROOT}/pointD.vertexCentred" \
            "${LEGACY_POINT_D_MAX}" "${LEGACY_POINT_D_MEAN}" \
            "${VERTEX_CENTRED_DISP_REL_TOL}"
        then
            failures=$((failures + 1))
        fi

        if ! run_vertex_centred_explicit_check; then
            failures=$((failures + 1))
        fi
    fi

    # The band above is a correctness bound. The same problem with the same
    # material constants must also give the legacy answer, far more closely
    # than the band allows
    if [[ -n "${UNSCOUPLED_EPS}" ]]
    then
        if awk "BEGIN {d = ${UNSCOUPLED_EPS} - ${LEGACY_UNSCOUPLED_EPS};
                       if (d < 0) d = -d;
                       exit !(d <= ${UNSCOUPLED_REL_TOL} * ${LEGACY_UNSCOUPLED_EPS})}"
        then
            printf "PASS: unsCoupled matches the legacy model (%.8g vs %.8g)\n" \
                "${UNSCOUPLED_EPS}" "${LEGACY_UNSCOUPLED_EPS}"
        else
            printf "FAIL: unsCoupled differs from the legacy model (%.8g vs %.8g)\n" \
                "${UNSCOUPLED_EPS}" "${LEGACY_UNSCOUPLED_EPS}"
            failures=$((failures + 1))
        fi
    fi
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
    # The sigma0 arms run the default petscSnes formulation, and the
    # cantileverTraction library they need is built by that arm's Allrun
    if [[ -z "${PETSC_DIR:-}" ]]; then
        echo "SKIP: sigma0 comparisons (PETSc is not installed)"
    elif make_sigma0_baseline; then
        for sigma0_mode in dict field both; do
            if ! run_sigma0_check "${sigma0_mode}"; then
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
