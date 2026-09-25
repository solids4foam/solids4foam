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
# crackingPlateHole regression test
#
# A plate with a hole, pulled at its left edge, cracks along its symmetry
# plane. simpleCrackerFvMesh does not change the mesh topology: the "down"
# patch is a symmetry plane whose faces the simpleCohesiveZone condition
# releases one at a time, under the Dugdale cohesive law, once their normal
# traction reaches sigmaMax. The quantities checked are the ones a cracking
# solve is about, none of them tied to a cell:
#
#   - the number of faces released, summed from the counts the mesh logs each
#     time step, and checked against the faces the condition logged releasing
#   - the time at which the reaction force on the loaded edge first falls,
#     which is where the crack first grows fast enough to shed load, and the
#     peak force reached just before it. This is the check most sensitive to
#     the cohesive law: 10 % off sigmaMax moves it by 2.7 %, where it moves the
#     final force by only 0.2 % and releases the same 30 faces
#   - the final reaction force on the loaded edge
#
# The case runs in seconds, so it is run to its full end time
# ============================================================

END_TIME=20

# Measured on foam-extend 4.1: 30 faces released, the force first falls at
# t = 8 from a peak force_y of 519120 at t = 7, final force_y 877672. The face
# count may move by two either way, the forces by 1 %
RELEASED_FACES_MIN=28
RELEASED_FACES_MAX=32
FIRST_DROP_TIME=8
PEAK_FORCE_MIN=5.14e5
PEAK_FORCE_MAX=5.24e5
FINAL_FORCE_MIN=8.69e5
FINAL_FORCE_MAX=8.86e5

# The answer of the removed legacy mechanicalModel, from the last commit that
# had it (mcl-stage8-coverage, c3a92b3d), on foam-extend 4.1: the faces
# released in each time step, as time:count, and the force_y history as time
# and force, written to twelve figures. The case is held to the history
# row by row, the largest difference relative to the legacy final force.
#
# The framework reproduced it exactly: the two force histories were identical
# to the 12 significant figures written, and the same faces were released at
# the same time steps. The two paths are the same linear elastic law, so any
# difference is round-off. Releasing a face is a threshold, so a real
# difference in the stress on the cohesive patch changes which face is
# released when, which moves the force by per cent. 1e-6 leaves room for
# round-off across compilers and is far below that
LEGACY_RELEASED_BY_STEP="3:3 4:3 5:3 6:2 7:2 8:4 9:2 10:3 11:1 12:2 13:1 14:1 15:1 16:1 17:1"
LEGACY_FORCE_HISTORY="
0 0
1 86934.8180971
2 173869.638267
3 258287.274213
4 333715.069706
5 401113.193395
6 462308.852288
7 519120.410696
8 498588.761597
9 509851.766338
10 523503.096582
11 549342.282304
12 576334.925567
13 612250.12715
14 648349.175807
15 684747.821768
16 721528.393483
17 758712.941656
18 796434.457649
19 834712.529485
20 877672.389608
"
LEGACY_FORCE_REL_TOL=1e-6

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
FORCE_FILE="postProcessing/0/solidForcesleft.dat"

echo "============================================================"
echo "crackingPlateHole regression test"
echo "Faces released in [${RELEASED_FACES_MIN}, ${RELEASED_FACES_MAX}]"
echo "Force_y first falls at t = ${FIRST_DROP_TIME}, from a peak in [${PEAK_FORCE_MIN}, ${PEAK_FORCE_MAX}]"
echo "Final force_y in [${FINAL_FORCE_MIN}, ${FINAL_FORCE_MAX}]"
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

    # Enough digits for the arm comparison to be about the solution rather
    # than about the last digit written
    sed -i 's/^writePrecision[[:space:]].*/writePrecision  12;/' "${controlDict}"

    # The reaction force on the loaded edge. The tutorial has no function
    # objects, so the copy gets one
    cat >> "${controlDict}" << 'EOF'

functions
{
    leftForce
    {
        type            solidForces;
        historyPatch    left;
    }
}
EOF
}

# The faces released at each time step, as "time<TAB>count" lines
released_per_step() {
    awk '/^Time = /{t = $3} /^Breaking [0-9]+ faces/{print t "\t" $2}' \
        "$1/${SOLVER_LOGFILE}"
}

# The faces released in each time step that released any, as "time:count"
# words on one line
released_by_step() {
    released_per_step "$1" | awk '
        {
            if (!($1 in s)) order[++n] = $1
            s[$1] += $2
        }
        END {
            for (i = 1; i <= n; i++) if (s[order[i]] > 0) printf "%s:%d ", order[i], s[order[i]]
            print ""
        }'
}

released_faces() {
    released_per_step "$1" | awk '{n += $2} END {print n + 0}'
}

# The faces the cohesive condition logged releasing
switched_faces() {
    grep -c "^Switching valueFraction to zero for face" \
        "$1/${SOLVER_LOGFILE}" || true
}

# The first time at which force_y is below its value at the previous time,
# and that previous value, as "time<TAB>force"
first_drop() {
    awk '!/^#/ && NF >= 3 {if (n++ && $3 < f) {print $1 "\t" f; exit} f = $3}' \
        "$1/${FORCE_FILE}"
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

# The case only runs on foam-extend, where simpleCrackerFvMesh is built
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
case_final_force=""

if ! arm_completed "${CASE_DIR}"; then
    fail "the case did not complete and converge"
else
    case_time=$(solids4Foam::latestTime "${CASE_DIR}")

    if [[ -z "${case_time}" ]] \
        || ! awk "BEGIN {exit !((${case_time} - ${END_TIME})^2 <= 1e-20)}"
    then
        fail "the case stopped at '${case_time}'; expected ${END_TIME}"
    fi

    case_faces=$(released_faces "${CASE_DIR}")
    switched=$(switched_faces "${CASE_DIR}")

    if (( case_faces != switched )); then
        fail "${case_faces} faces released but ${switched} logged by the cohesive condition"
    elif in_range "${case_faces}" \
        "${RELEASED_FACES_MIN}" "${RELEASED_FACES_MAX}"
    then
        echo "PASS: faces released = ${case_faces}"
    else
        fail "faces released = ${case_faces}"
    fi

    if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
        fail "the case wrote no ${FORCE_FILE}"
    else
        drop_time=""
        peak_force=""
        read -r drop_time peak_force <<< "$(first_drop "${CASE_DIR}")" || true
        case_final_force=$(final_force "${CASE_DIR}")

        if [[ -n "${drop_time}" && -n "${peak_force}" ]] \
         && awk "BEGIN {exit !((${drop_time} - ${FIRST_DROP_TIME})^2 <= 1e-20)}" \
         && in_range "${peak_force}" "${PEAK_FORCE_MIN}" "${PEAK_FORCE_MAX}"
        then
            printf "PASS: force_y first falls at t = %s, from %.6g\n" \
                "${drop_time}" "${peak_force}"
        else
            fail "force_y first falls at t = '${drop_time}', from '${peak_force}'"
        fi

        if [[ -n "${case_final_force}" ]] \
         && in_range "${case_final_force}" \
                "${FINAL_FORCE_MIN}" "${FINAL_FORCE_MAX}"
        then
            printf "PASS: final force_y = %.6g\n" "${case_final_force}"
        else
            fail "final force_y = '${case_final_force}'"
        fi
    fi
fi

# ------------------------------------------------------------
# Against the legacy answer
#
# simpleCohesiveZone sets its traction through the solid model's
# tractionBoundarySnGrad, which takes the implicit stiffness from impK, so
# this runs on lawImpK()
# ------------------------------------------------------------

if grep -q "Selecting mechanical constitutive law" \
    "${CASE_DIR}/${SOLVER_LOGFILE}" 2>/dev/null
then
    echo "PASS: the material came from the framework"
else
    fail "the case constructed no mechanical constitutive law"
fi

if [[ -n "${case_time}" && -f "${CASE_DIR}/${FORCE_FILE}" ]]; then
    # Releasing is a threshold, so the same answer releases the same faces at
    # the same time steps
    released="$(released_by_step "${CASE_DIR}")"
    if [[ "${released% }" == "${LEGACY_RELEASED_BY_STEP}" ]]
    then
        echo "PASS: the same faces released at the same times as the legacy model"
    else
        fail "faces released at different times from the legacy model ('${released% }')"
    fi

    # Row by row, over the whole history, the case having written the times
    # the legacy model did, relative to the legacy final force
    if rel=$(awk '
            FNR == NR {
                if (NF == 2) {f[$1] = $2; n++; ref = $2}
                next
            }
            !/^#/ && NF >= 3 {
                if (!($1 in f)) {bad = 1; exit}
                d = $3 - f[$1]; if (d < 0) d = -d
                if (d > m) m = d
                k++
            }
            END {if (bad || k != n || n == 0) exit 1; printf "%.6g", m/ref}
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
