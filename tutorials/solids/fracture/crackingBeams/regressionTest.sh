#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# crackingBeams regression test
#
# Two beams are pulled apart at their left ends and crack along their
# mid-planes. crackerFvMesh breaks internal faces onto the "crack" patch, whose
# faces carry the variableMixedMode cohesive zone law. The quantities checked
# are the ones a cracking solve is about, none of them tied to a cell, since
# the mesh topology changes as it cracks:
#
#   - the number of faces on the crack patch at the end, read from the mesh
#     the solver writes, and checked against the number of face breaks it
#     logged: each break puts two faces on the patch, one per crack flank
#   - the peak reaction force on the loaded patch, and the time it is reached:
#     the force rises elastically, peaks as the cohesive zone first gives way,
#     then falls as the crack runs
#   - the final reaction force, well below the peak, which is only true if the
#     crack has propagated and the beams have softened
#
# The run is shortened from 20 to 10 time steps: the peak is at step 6, and by
# step 10 the force has fallen by a quarter, so every check above still has
# something to measure, in half the time
# ============================================================

END_TIME=10

# Measured on foam-extend 4.1: 44 crack faces (22 breaks), peak force_y
# 81.6217 N at t = 6, final force_y 61.4155 N. The face count may move by one
# break either way, the forces by 1 %
CRACK_FACES_MIN=40
CRACK_FACES_MAX=48
PEAK_TIME=6
PEAK_FORCE_MIN=80.8
PEAK_FORCE_MAX=82.4
FINAL_FORCE_MIN=60.8
FINAL_FORCE_MAX=62.0

# The answer of the removed legacy mechanicalModel, from the last commit that
# had it (mcl-stage8-coverage, c3a92b3d), on foam-extend 4.1: the crack patch
# size, and the force_y history as time and force, written to twelve figures.
# The case is held to the history row by row, the largest difference relative
# to the legacy peak.
#
# The framework reproduced it exactly: the two force histories were identical
# to the 12 significant figures written, with the same iteration counts and
# the same faces broken at the same iterations. The two paths are the same
# linear elastic law, so any difference is round-off. The bound matters more
# than it would on a smooth problem: breaking a face is a threshold, and a
# difference in the implicit stiffness the cohesive law takes its penalty
# from, or in the stress on the new crack faces, changes which face breaks
# when, which moves the force by per cent. 1e-6 leaves room for round-off
# across compilers and is far below that
LEGACY_CRACK_FACES=44
LEGACY_FORCE_HISTORY="
0 0
1 26.1319882895
2 44.4955884866
3 58.1396766671
4 69.2483734867
5 78.671274439
6 81.6216548573
7 73.9834926861
8 71.0552239168
9 64.2618226861
10 61.4154571245
"
LEGACY_FORCE_REL_TOL=1e-6

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
CRACK_PATCH="crack"
FORCE_FILE="postProcessing/0/solidForcestopLoading.dat"

echo "============================================================"
echo "crackingBeams regression test"
echo "Crack patch faces at t = ${END_TIME} in [${CRACK_FACES_MIN}, ${CRACK_FACES_MAX}]"
echo "Peak force_y at t = ${PEAK_TIME} in [${PEAK_FORCE_MIN}, ${PEAK_FORCE_MAX}] N"
echo "Final force_y in [${FINAL_FORCE_MIN}, ${FINAL_FORCE_MAX}] N"
echo "force_y history against the legacy model, rel. diff <= ${LEGACY_FORCE_REL_TOL}"
echo "============================================================"
echo

prepare_arm() {
    local d="$1"

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

    local controlDict="${d}/system/controlDict"

    sed -i 's/^endTime[[:space:]]\+20;/endTime         '"${END_TIME}"';/' \
        "${controlDict}"

    # Enough digits for the arm comparison to be about the solution rather
    # than about the last digit written
    sed -i 's/^writePrecision[[:space:]].*/writePrecision  12;/' "${controlDict}"

    # The reaction force on the loaded patch. The tutorial has no function
    # objects, so the copy gets one
    cat >> "${controlDict}" << 'EOF'

functions
{
    topLoadingForce
    {
        type            solidForces;
        historyPatch    topLoading;
    }
}
EOF
}

# The number of faces on the crack patch in the latest mesh the case wrote
crack_patch_faces() {
    local d="$1"
    local t
    t=$(solids4Foam::latestTime "${d}")
    [[ -n "${t}" && -f "${d}/${t}/polyMesh/boundary" ]] || return 0

    awk -v patch="${CRACK_PATCH}" '
        $1 == patch {inPatch = 1; next}
        inPatch && $1 == "nFaces" {sub(";", "", $2); print $2; exit}
        inPatch && $1 == "}" {inPatch = 0}
    ' "${d}/${t}/polyMesh/boundary"
}

# The number of internal faces the solver logged breaking
logged_breaks() {
    grep -c "^Breaking internal face" "$1/${SOLVER_LOGFILE}" || true
}

# The peak force_y and the time it is reached, as "time<TAB>force"
peak_force() {
    awk '!/^#/ && NF >= 3 {if (!n++ || $3 > f) {f = $3; t = $1}}
         END {if (n) print t "\t" f}' "$1/${FORCE_FILE}"
}

final_force() {
    awk '!/^#/ && NF >= 3 {f = $3} END {print f}' "$1/${FORCE_FILE}"
}

in_range() {
    awk "BEGIN {exit !($1 >= $2 && $1 <= $3)}"
}

arm_completed() {
    local d="$1"
    grep -q "^End" "${d}/${SOLVER_LOGFILE}" 2>/dev/null \
        && ! grep -qE "FOAM FATAL|Maximum iterations reached" \
            "${d}/${SOLVER_LOGFILE}"
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
    prepare_arm "${CASE_DIR}"

    # A failed run is reported by the checks below rather than by set -e, so
    # that the case is still cleaned up
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

clean_arms() {
    if [ "$CHECK_ONLY" = false ]; then
        ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    fi
}

# The case only runs on foam-extend, where crackerFvMesh is built
if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    clean_arms
    exit 0
fi

failures=0

fail() {
    echo "FAIL: $*"
    failures=$((failures + 1))
}

# ------------------------------------------------------------
# Against the measured answer
# ------------------------------------------------------------

case_time=""
case_faces=""
case_peak_force=""

if ! arm_completed "${CASE_DIR}"; then
    fail "the case did not complete and converge"
else
    case_time=$(solids4Foam::latestTime "${CASE_DIR}")

    if [[ -z "${case_time}" ]] \
        || ! awk "BEGIN {exit !((${case_time} - ${END_TIME})^2 <= 1e-20)}"
    then
        fail "the case stopped at '${case_time}'; expected ${END_TIME}"
    fi

    case_faces=$(crack_patch_faces "${CASE_DIR}")
    breaks=$(logged_breaks "${CASE_DIR}")

    if [[ -z "${case_faces}" ]]; then
        fail "could not read the ${CRACK_PATCH} patch size"
    elif (( case_faces != 2*breaks )); then
        fail "${CRACK_PATCH} patch has ${case_faces} faces but ${breaks} breaks were logged"
    elif in_range "${case_faces}" "${CRACK_FACES_MIN}" "${CRACK_FACES_MAX}"
    then
        echo "PASS: ${CRACK_PATCH} patch faces = ${case_faces} (${breaks} breaks)"
    else
        fail "${CRACK_PATCH} patch faces = ${case_faces} (${breaks} breaks)"
    fi

    if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
        fail "the case wrote no ${FORCE_FILE}"
    else
        read -r peak_time case_peak_force <<< "$(peak_force "${CASE_DIR}")"
        case_final_force=$(final_force "${CASE_DIR}")

        if [[ -z "${case_peak_force}" || -z "${case_final_force}" ]]; then
            fail "could not extract the force history"
        else
            if in_range "${case_peak_force}" \
                "${PEAK_FORCE_MIN}" "${PEAK_FORCE_MAX}" \
             && awk "BEGIN {exit !((${peak_time} - ${PEAK_TIME})^2 <= 1e-20)}"
            then
                printf "PASS: peak force_y = %.6g N at t = %s\n" \
                    "${case_peak_force}" "${peak_time}"
            else
                fail "peak force_y = ${case_peak_force} N at t = ${peak_time}"
            fi

            if in_range "${case_final_force}" \
                "${FINAL_FORCE_MIN}" "${FINAL_FORCE_MAX}"
            then
                printf "PASS: final force_y = %.6g N\n" "${case_final_force}"
            else
                fail "final force_y = ${case_final_force} N"
            fi
        fi
    fi
fi

# ------------------------------------------------------------
# Against the legacy answer
#
# The cohesive zone laws take their penalty stiffness from the field the solid
# model registers as impK, looked up by name, which is lawImpK(), so this
# is the test of it for fracture, as curvedBeams is for contact and 3dTube is
# for fluid-solid interaction. It is also the framework's run on a mesh whose
# topology changes: the crack patch gains faces as it cracks, and the
# framework's boundary addressing and the solid model's density have to follow
# ------------------------------------------------------------

if grep -q "Selecting mechanical constitutive law" \
    "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null
then
    echo "PASS: the material came from the framework"
else
    fail "the case constructed no mechanical constitutive law"
fi

if [[ -n "${case_faces}" && -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
    # Breaking is a threshold, so the same answer breaks exactly the same
    # faces
    if [[ "${case_faces}" == "${LEGACY_CRACK_FACES}" ]]; then
        echo "PASS: ${CRACK_PATCH} patch faces = ${case_faces}, as the legacy model"
    else
        fail "${CRACK_PATCH} patch faces = '${case_faces}', legacy model ${LEGACY_CRACK_FACES}"
    fi

    # Row by row, over the whole history, the case having written the times
    # the legacy model did
    if rel=$(awk '
            FNR == NR {
                if (NF == 2) {f[$1] = $2; n++; if ($2 > peak) peak = $2}
                next
            }
            !/^#/ && NF >= 3 {
                if (!($1 in f)) {bad = 1; exit}
                d = $3 - f[$1]; if (d < 0) d = -d
                if (d > m) m = d
                k++
            }
            END {if (bad || k != n || n == 0) exit 1; printf "%.6g", m/peak}
        ' <(echo "${LEGACY_FORCE_HISTORY}") "${CASE_DIR}/${FORCE_FILE}")
    then
        if awk "BEGIN {exit !(${rel} <= ${LEGACY_FORCE_REL_TOL})}"; then
            echo "PASS: force_y history matches the legacy model, max rel. diff = ${rel}"
        else
            fail "force_y history differs from the legacy model, max rel. diff = ${rel}"
        fi
    else
        fail "the force history has different times from the legacy model's"
    fi
fi

clean_arms

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
