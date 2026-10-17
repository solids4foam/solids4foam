#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
source "${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

# Check one full oscillation for BDF2, trapezoidal Newmark and Bossak-Newmark.
# The coarse-mesh peak is close to 2.725 m for all three schemes on v2512.
# The BDF2 peaks on foam-extend-4.1 and OpenFOAM-9 are 2.7404 and 2.7556 m.
# Keep the existing cross-version amplitude band; the operator check also
# verifies the selected scheme and its physical/weighted accelerations.
PEAK_MIN=2.70
PEAK_MAX=2.78
END_TIME=0.65
SCHEMES=(bdf2 newmark bossak)
CHECK_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run) CHECK_ONLY=true ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

prepare_case() {
    local case_dir="$1"
    rm -rf "$case_dir"
    mkdir -p "$case_dir"
    cp -a "${SCRIPT_DIR}/0" "${SCRIPT_DIR}/constant" \
        "${SCRIPT_DIR}/system" "${SCRIPT_DIR}/reference" \
        "${SCRIPT_DIR}/Allrun" "${SCRIPT_DIR}/Allclean" \
        "${SCRIPT_DIR}/plot.gnuplot" "$case_dir/"
}

failures=0
for scheme in "${SCHEMES[@]}"; do
    case_dir="${REGRESSION_ROOT}/${scheme}"
    echo "============================================================"
    echo "cantileverVibration: ${scheme} (PETSc SNES)"
    echo "Peak tip displacement magnitude in [${PEAK_MIN}, ${PEAK_MAX}] m"

    if [ "$CHECK_ONLY" = false ]; then
        prepare_case "$case_dir"
        if ! (cd "$case_dir" && ./Allclean > /dev/null 2>&1 \
            && ./Allrun petscSnes "$scheme" > log.Allrun 2>&1)
        then
            echo "FAIL: Allrun failed; see ${case_dir}/log.Allrun"
            failures=$((failures + 1))
            continue
        fi
    fi

    if solids4Foam::regressionCaseSkipped "${case_dir}/log.Allrun"; then
        echo "Skipping regression checks: tutorial skipped in this environment"
        exit 0
    fi

    solver_log="${case_dir}/log.solids4Foam"
    if [[ ! -f "$solver_log" ]] \
        || grep -Eq 'FOAM FATAL|^ERROR$|\[stack trace\]' "$solver_log" \
        || ! grep -q '^End' "$solver_log"
    then
        echo "FAIL: solids4Foam did not complete; see ${solver_log}"
        failures=$((failures + 1))
        continue
    fi

    disp_file="${case_dir}/postProcessing/0/solidPointDisplacement_pointDisp.dat"
    if [[ ! -s "$disp_file" ]]; then
        echo "FAIL: Missing point displacement output"
        failures=$((failures + 1))
        continue
    fi

    if ! awk -v end="$END_TIME" '
        !/^#/ && NF >= 5 {t = $1; n++}
        END {exit !(n && t >= end - 1e-10)}' "$disp_file"
    then
        echo "FAIL: Displacement history does not reach ${END_TIME} s"
        failures=$((failures + 1))
        continue
    fi

    peak=$(awk '!/^#/ && NF >= 5 {if (!n++ || $5 > m) m = $5}
        END {if (n) print m}' "$disp_file")
    if awk -v peak="$peak" -v low="$PEAK_MIN" -v high="$PEAK_MAX" \
        'BEGIN {exit !(peak >= low && peak <= high)}'
    then
        echo "PASS: Peak tip displacement = ${peak} m"
    else
        echo "FAIL: Peak tip displacement = ${peak} m"
        failures=$((failures + 1))
    fi

    # Use the same mesh and scheme dictionaries for the operator-level test.
    # Its manufactured displacement does not read or change the solver fields.
    if [ "$CHECK_ONLY" = false ]; then
        if ! (cd "$case_dir" && Test-fvcD2dt2 > log.Test-fvcD2dt2 2>&1); then
            echo "FAIL: Test-fvcD2dt2; see ${case_dir}/log.Test-fvcD2dt2"
            failures=$((failures + 1))
            continue
        fi
    fi
    case "$scheme" in
        bdf2) selected='Selected d2dt2 scheme: backward' ;;
        newmark) selected='Selected d2dt2 scheme: NewmarkBeta beta=0.25 gamma=0.5 alphaM=0' ;;
        bossak) selected='Selected d2dt2 scheme: NewmarkBeta beta=0.3025 gamma=0.6 alphaM=-0.1' ;;
    esac
    if grep -Fxq "$selected" "${case_dir}/log.Test-fvcD2dt2" \
        && grep -q '^Test-fvcD2dt2: PASSED' "${case_dir}/log.Test-fvcD2dt2"
    then
        echo "PASS: ${scheme} inertia and physical-acceleration checks"
    else
        echo "FAIL: Missing operator-test success or incorrect time scheme"
        failures=$((failures + 1))
    fi

    # newmark selects NewmarkBeta in d2dt2Schemes only, so the linear
    # predictor reports that it uses the ddtSchemes default; bossak also has
    # the optional matching ddtSchemes entry, so it does not
    predictor_note='uses the ddtSchemes scheme backward'
    case "$scheme" in
        newmark)
            if grep -Fq "$predictor_note" "$solver_log"; then
                echo "PASS: d2dt2Schemes-only NewmarkBeta accepted"
            else
                echo "FAIL: Missing the ddtSchemes predictor note"
                failures=$((failures + 1))
            fi
            ;;
        bossak)
            if ! grep -Fq "$predictor_note" "$solver_log"; then
                echo "PASS: Matching ddtSchemes NewmarkBeta entry accepted"
            else
                echo "FAIL: The matching ddtSchemes entry was not used"
                failures=$((failures + 1))
            fi
            ;;
    esac
done

# A NewmarkBeta ddtSchemes entry whose coefficients differ from d2dt2Schemes
# must stop with a fatal error, as the two would advance the same stored
# state. Test-fvcD2dt2 calls fvm::d2dt2, which makes the check, on the bossak
# mesh with trapezoidal coefficients in ddtSchemes
mismatch_dir="${REGRESSION_ROOT}/mismatch"
mismatch_log="${mismatch_dir}/log.Test-fvcD2dt2"
echo "============================================================"
echo "cantileverVibration: mismatched NewmarkBeta coefficients"

if [ "$CHECK_ONLY" = false ]; then
    rm -rf "$mismatch_dir"
    mkdir -p "$mismatch_dir/system"
    cp -a "${REGRESSION_ROOT}/bossak/constant" "$mismatch_dir/"
    cp -a "${SCRIPT_DIR}/system/controlDict" \
        "${SCRIPT_DIR}/system/fvSolution.petscSnes" "$mismatch_dir/system/"
    mv "$mismatch_dir/system/fvSolution.petscSnes" \
        "$mismatch_dir/system/fvSolution"
    awk '
        /^ddtSchemes/ {inDdt = 1}
        inDdt && /NewmarkBeta/ {sub(/NewmarkBeta.*;/, "NewmarkBeta 0.25 0.5;")}
        inDdt && /^}/ {inDdt = 0}
        {print}' "${SCRIPT_DIR}/system/fvSchemes.bossak" \
        > "$mismatch_dir/system/fvSchemes"

    if (cd "$mismatch_dir" && Test-fvcD2dt2 > log.Test-fvcD2dt2 2>&1); then
        echo "FAIL: Test-fvcD2dt2 ran with mismatched coefficients"
        failures=$((failures + 1))
    fi
fi

if [[ -f "$mismatch_log" ]] \
    && grep -q 'selects it with different coefficients' "$mismatch_log"
then
    echo "PASS: Mismatched coefficients stop with a fatal error"
else
    echo "FAIL: No mismatched-coefficient error; see ${mismatch_log}"
    failures=$((failures + 1))
fi

if (( failures == 0 )); then
    echo "Regression test PASSED"
else
    echo "Regression test FAILED (${failures} checks)"
    exit 1
fi
