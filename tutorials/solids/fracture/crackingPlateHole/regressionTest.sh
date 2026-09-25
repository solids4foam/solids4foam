#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
LEGACY_DIR="${REGRESSION_ROOT}/legacy"
FRAMEWORK_DIR="${REGRESSION_ROOT}/framework"
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

# The framework arm against the legacy arm: the largest difference in force_y
# over the whole history, relative to the legacy final force.
#
# Measured at zero: the two force histories are identical to the 12
# significant figures written, and the same faces are released at the same
# time steps. The two paths are the same linear elastic law, so any difference
# is round-off. Releasing a face is a threshold, so a real difference in the
# stress on the cohesive patch changes which face is released when, which
# moves the force by per cent. 1e-6 leaves room for round-off across compilers
# and is far below that
FRAMEWORK_FORCE_REL_TOL=1e-6

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
FORCE_FILE="postProcessing/0/solidForcesleft.dat"

echo "============================================================"
echo "crackingPlateHole regression test"
echo "Faces released in [${RELEASED_FACES_MIN}, ${RELEASED_FACES_MAX}]"
echo "Force_y first falls at t = ${FIRST_DROP_TIME}, from a peak in [${PEAK_FORCE_MIN}, ${PEAK_FORCE_MAX}]"
echo "Final force_y in [${FINAL_FORCE_MIN}, ${FINAL_FORCE_MAX}]"
echo "Framework vs legacy force_y rel. diff <= ${FRAMEWORK_FORCE_REL_TOL}"
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
    prepare_arm "${LEGACY_DIR}"
    prepare_arm "${FRAMEWORK_DIR}"

    # The two arms differ in this one entry and nothing else. It goes inside
    # the solid model's coeffs sub-dictionary, which is where the model reads
    # it
    sed -i \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${FRAMEWORK_DIR}/constant/solidProperties"

    if ! grep -q "useMechanicalConstitutiveLawManager" \
        "${FRAMEWORK_DIR}/constant/solidProperties"
    then
        echo "FAIL: could not set the framework switch on the framework arm"
        exit 1
    fi

    # A failed run is reported by the checks below rather than by set -e, so
    # that both arms are run and cleaned up
    ( cd "${LEGACY_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true
    ( cd "${FRAMEWORK_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

clean_arms() {
    if [ "$CHECK_ONLY" = false ]; then
        local d
        for d in "${LEGACY_DIR}" "${FRAMEWORK_DIR}"; do
            ( cd "${d}" && ./Allclean > /dev/null 2>&1 ) || true
        done
    fi
}

# The case only runs on foam-extend, where simpleCrackerFvMesh is built
if solids4Foam::regressionCaseSkipped "${LEGACY_DIR}/${ALLRUN_LOGFILE}"; then
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
# The legacy arm against the measured answer
# ------------------------------------------------------------

legacy_time=""
legacy_faces=""
legacy_final_force=""

if ! arm_completed "${LEGACY_DIR}"; then
    fail "the legacy arm did not complete and converge"
else
    legacy_time=$(solids4Foam::latestTime "${LEGACY_DIR}")

    if [[ -z "${legacy_time}" ]] \
        || ! awk "BEGIN {exit !((${legacy_time} - ${END_TIME})^2 <= 1e-20)}"
    then
        fail "the legacy arm stopped at '${legacy_time}'; expected ${END_TIME}"
    fi

    legacy_faces=$(released_faces "${LEGACY_DIR}")
    switched=$(switched_faces "${LEGACY_DIR}")

    if (( legacy_faces != switched )); then
        fail "${legacy_faces} faces released but ${switched} logged by the cohesive condition"
    elif in_range "${legacy_faces}" \
        "${RELEASED_FACES_MIN}" "${RELEASED_FACES_MAX}"
    then
        echo "PASS: faces released = ${legacy_faces}"
    else
        fail "faces released = ${legacy_faces}"
    fi

    if [[ ! -f "${LEGACY_DIR}/${FORCE_FILE}" ]]; then
        fail "the legacy arm wrote no ${FORCE_FILE}"
    else
        drop_time=""
        peak_force=""
        read -r drop_time peak_force <<< "$(first_drop "${LEGACY_DIR}")" || true
        legacy_final_force=$(final_force "${LEGACY_DIR}")

        if [[ -n "${drop_time}" && -n "${peak_force}" ]] \
         && awk "BEGIN {exit !((${drop_time} - ${FIRST_DROP_TIME})^2 <= 1e-20)}" \
         && in_range "${peak_force}" "${PEAK_FORCE_MIN}" "${PEAK_FORCE_MAX}"
        then
            printf "PASS: force_y first falls at t = %s, from %.6g\n" \
                "${drop_time}" "${peak_force}"
        else
            fail "force_y first falls at t = '${drop_time}', from '${peak_force}'"
        fi

        if [[ -n "${legacy_final_force}" ]] \
         && in_range "${legacy_final_force}" \
                "${FINAL_FORCE_MIN}" "${FINAL_FORCE_MAX}"
        then
            printf "PASS: final force_y = %.6g\n" "${legacy_final_force}"
        else
            fail "final force_y = '${legacy_final_force}'"
        fi
    fi
fi

# ------------------------------------------------------------
# The framework arm against the legacy arm
#
# simpleCohesiveZone sets its traction through the solid model's
# tractionBoundarySnGrad, which takes the implicit stiffness from impK, so on
# the framework arm this runs on frameworkImpK()
# ------------------------------------------------------------

# Each arm must have taken the path it was set up for, or the comparison is a
# run against itself
if grep -q "Selecting mechanical constitutive law" \
    "${FRAMEWORK_DIR}/${SOLVER_LOGFILE}" 2>/dev/null
then
    echo "PASS: framework arm took the framework path"
else
    fail "the framework arm did not use the framework"
fi

if grep -q "Selecting mechanical constitutive law" \
    "${LEGACY_DIR}/${SOLVER_LOGFILE}" 2>/dev/null
then
    fail "the legacy arm used the framework"
else
    echo "PASS: legacy arm took the legacy path"
fi

if ! arm_completed "${FRAMEWORK_DIR}"; then
    fail "the framework arm did not complete and converge"
elif [[ -n "${legacy_time}" && -n "${legacy_final_force}" ]]; then
    framework_time=$(solids4Foam::latestTime "${FRAMEWORK_DIR}")

    if [[ "${framework_time}" != "${legacy_time}" ]]; then
        fail "the arms reached different times ('${legacy_time}' vs '${framework_time}')"
    fi

    # Releasing is a threshold, so the same answer releases the same faces at
    # the same time steps
    if [[ "$(released_per_step "${FRAMEWORK_DIR}")" \
        == "$(released_per_step "${LEGACY_DIR}")" ]]
    then
        echo "PASS: framework arm released the same faces at the same times"
    else
        fail "framework arm released faces at different times ($(released_faces "${FRAMEWORK_DIR}") faces, legacy ${legacy_faces})"
    fi

    if [[ ! -f "${FRAMEWORK_DIR}/${FORCE_FILE}" ]]; then
        fail "the framework arm wrote no ${FORCE_FILE}"
    else
        # Row by row, over the whole history, both arms having written the
        # same times
        if rel=$(awk -v ref="${legacy_final_force}" '
                !/^#/ && NF >= 3 {
                    if (FNR == NR) {f[$1] = $3; n++; next}
                    if (!($1 in f)) {bad = 1; exit}
                    d = $3 - f[$1]; if (d < 0) d = -d
                    if (d > m) m = d
                    k++
                }
                END {if (bad || k != n || n == 0) exit 1; printf "%.6g", m/ref}
            ' "${LEGACY_DIR}/${FORCE_FILE}" "${FRAMEWORK_DIR}/${FORCE_FILE}")
        then
            if awk "BEGIN {exit !(${rel} <= ${FRAMEWORK_FORCE_REL_TOL})}"; then
                echo "PASS: framework and legacy force_y agree, max rel. diff = ${rel}"
            else
                fail "framework and legacy force_y differ, max rel. diff = ${rel}"
            fi
        else
            fail "the two arms' force histories have different times"
        fi
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
