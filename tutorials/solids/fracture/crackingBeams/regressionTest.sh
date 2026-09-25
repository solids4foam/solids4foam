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

# The framework arm against the legacy arm: the largest difference in force_y
# over the whole history, relative to the legacy peak.
#
# Measured at zero: the two force histories are identical to the 12
# significant figures written, with the same iteration counts and the same
# faces broken at the same iterations. The two paths are the same linear
# elastic law, so any difference is round-off. The bound matters more
# than it would on a smooth problem: breaking a face is a threshold, and a
# difference in the implicit stiffness the cohesive law takes its penalty from,
# or in the stress on the new crack faces, changes which face breaks when,
# which moves the force by per cent. 1e-6 leaves room for round-off across
# compilers and is far below that
FRAMEWORK_FORCE_REL_TOL=1e-6

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
CRACK_PATCH="crack"
FORCE_FILE="postProcessing/0/solidForcestopLoading.dat"

echo "============================================================"
echo "crackingBeams regression test"
echo "Crack patch faces at t = ${END_TIME} in [${CRACK_FACES_MIN}, ${CRACK_FACES_MAX}]"
echo "Peak force_y at t = ${PEAK_TIME} in [${PEAK_FORCE_MIN}, ${PEAK_FORCE_MAX}] N"
echo "Final force_y in [${FINAL_FORCE_MIN}, ${FINAL_FORCE_MAX}] N"
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

# The case only runs on foam-extend, where crackerFvMesh is built
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
legacy_peak_force=""

if ! arm_completed "${LEGACY_DIR}"; then
    fail "the legacy arm did not complete and converge"
else
    legacy_time=$(solids4Foam::latestTime "${LEGACY_DIR}")

    if [[ -z "${legacy_time}" ]] \
        || ! awk "BEGIN {exit !((${legacy_time} - ${END_TIME})^2 <= 1e-20)}"
    then
        fail "the legacy arm stopped at '${legacy_time}'; expected ${END_TIME}"
    fi

    legacy_faces=$(crack_patch_faces "${LEGACY_DIR}")
    breaks=$(logged_breaks "${LEGACY_DIR}")

    if [[ -z "${legacy_faces}" ]]; then
        fail "could not read the ${CRACK_PATCH} patch size"
    elif (( legacy_faces != 2*breaks )); then
        fail "${CRACK_PATCH} patch has ${legacy_faces} faces but ${breaks} breaks were logged"
    elif in_range "${legacy_faces}" "${CRACK_FACES_MIN}" "${CRACK_FACES_MAX}"
    then
        echo "PASS: ${CRACK_PATCH} patch faces = ${legacy_faces} (${breaks} breaks)"
    else
        fail "${CRACK_PATCH} patch faces = ${legacy_faces} (${breaks} breaks)"
    fi

    if [[ ! -f "${LEGACY_DIR}/${FORCE_FILE}" ]]; then
        fail "the legacy arm wrote no ${FORCE_FILE}"
    else
        read -r peak_time legacy_peak_force <<< "$(peak_force "${LEGACY_DIR}")"
        legacy_final_force=$(final_force "${LEGACY_DIR}")

        if [[ -z "${legacy_peak_force}" || -z "${legacy_final_force}" ]]; then
            fail "could not extract the force history"
        else
            if in_range "${legacy_peak_force}" \
                "${PEAK_FORCE_MIN}" "${PEAK_FORCE_MAX}" \
             && awk "BEGIN {exit !((${peak_time} - ${PEAK_TIME})^2 <= 1e-20)}"
            then
                printf "PASS: peak force_y = %.6g N at t = %s\n" \
                    "${legacy_peak_force}" "${peak_time}"
            else
                fail "peak force_y = ${legacy_peak_force} N at t = ${peak_time}"
            fi

            if in_range "${legacy_final_force}" \
                "${FINAL_FORCE_MIN}" "${FINAL_FORCE_MAX}"
            then
                printf "PASS: final force_y = %.6g N\n" "${legacy_final_force}"
            else
                fail "final force_y = ${legacy_final_force} N"
            fi
        fi
    fi
fi

# ------------------------------------------------------------
# The framework arm against the legacy arm
#
# The cohesive zone laws take their penalty stiffness from the field the solid
# model registers as impK, looked up by name. On the framework arm that is
# frameworkImpK(), so this is the test of it for fracture, as curvedBeams is
# for contact and 3dTube is for fluid-solid interaction. It is also the
# framework's first run on a mesh whose topology changes: the crack patch
# gains faces as it cracks, and the framework's boundary addressing and the
# solid model's density have to follow
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
elif [[ -n "${legacy_time}" && -n "${legacy_faces}" \
     && -n "${legacy_peak_force}" ]]
then
    framework_time=$(solids4Foam::latestTime "${FRAMEWORK_DIR}")
    framework_faces=$(crack_patch_faces "${FRAMEWORK_DIR}")

    if [[ "${framework_time}" != "${legacy_time}" ]]; then
        fail "the arms reached different times ('${legacy_time}' vs '${framework_time}')"
    fi

    # Breaking is a threshold, so the same answer breaks exactly the same
    # faces
    if [[ "${framework_faces}" == "${legacy_faces}" ]]; then
        echo "PASS: framework arm ${CRACK_PATCH} patch faces = ${framework_faces}"
    else
        fail "framework arm ${CRACK_PATCH} patch faces = '${framework_faces}', legacy ${legacy_faces}"
    fi

    if [[ ! -f "${FRAMEWORK_DIR}/${FORCE_FILE}" ]]; then
        fail "the framework arm wrote no ${FORCE_FILE}"
    else
        # Row by row, over the whole history, both arms having written the
        # same times
        if rel=$(awk -v peak="${legacy_peak_force}" '
                !/^#/ && NF >= 3 {
                    if (FNR == NR) {f[$1] = $3; n++; next}
                    if (!($1 in f)) {bad = 1; exit}
                    d = $3 - f[$1]; if (d < 0) d = -d
                    if (d > m) m = d
                    k++
                }
                END {if (bad || k != n || n == 0) exit 1; printf "%.6g", m/peak}
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
