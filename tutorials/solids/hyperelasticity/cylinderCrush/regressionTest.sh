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
# cylinderCrush regression test
# Uses the short displacement and force histories as a contact benchmark.
# ============================================================

FORCE_Y_MIN=-550
FORCE_Y_MAX=-545
DISP_Y_MIN=-0.0034
DISP_Y_MAX=-0.0032

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.solids4Foam"
FORCE_FILE="postProcessing/0/solidForcescylinderContact.dat"
DISP_FILE="postProcessing/0/solidPointDisplacement_displacement.dat"

echo "============================================================"
echo "cylinderCrush regression test"
echo "Final cylinder force_y in [${FORCE_Y_MIN}, ${FORCE_Y_MAX}] N"
echo "Final probe disp_y in [${DISP_Y_MIN}, ${DISP_Y_MAX}] m"
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

    sed -i.bak 's/^endTime[[:space:]]\+30;/endTime         1;/' "${CASE_DIR}/system/controlDict"
    rm -f "${CASE_DIR}/system/controlDict.bak"
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
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

if [[ ! -f "${CASE_DIR}/${FORCE_FILE}" || ! -f "${CASE_DIR}/${DISP_FILE}" ]]; then
    echo "FAIL: Could not find one or more history files"
    exit 1
fi

final_force_y=$(awk 'END {print $3}' "${CASE_DIR}/${FORCE_FILE}")
final_disp_y=$(awk 'END {print $3}' "${CASE_DIR}/${DISP_FILE}")

if [[ -z "${final_force_y}" || -z "${final_disp_y}" ]]; then
    echo "FAIL: Could not extract final force/displacement"
    exit 1
fi

# Run the case both ways and require the two to agree.
#
# Both arms turn the pressure equation off. The case ships with
# solvePressureEqn yes and a smoothing scale factor, and the framework has no
# equivalent yet - the hydrostatic stress there is 0.5*K*(J^2 - 1) taken
# directly, where the legacy law solves and smooths an equation for it. So this
# checks the Ogden law itself: the principal stretches, the stress built from
# them, and the rotation back. The shipped configuration is covered by the
# checks above, and the pressure equation is a separate piece of work.
#
# Two of the thirty steps. The legacy solver reaches its corrector limit from
# the third step onwards, so past that it is an unconverged answer being
# compared against a converged one
COMPARISON_END_TIME=2

# The latest written time directory.
#
# Not foamListTimes: it needs an etc/controlDict that this foam-extend
# installation does not provide, and this case runs only on foam-extend
latest_time_dir() {
    ls -1 "$1" 2>/dev/null \
        | grep -E '^[0-9]+([.][0-9]+)?$' \
        | sort -g \
        | tail -n 1
}

run_framework_comparison() {
    local legacy_dir="${REGRESSION_ROOT}/comparisonLegacy"
    local framework_dir="${REGRESSION_ROOT}/comparisonFramework"
    local dir

    for dir in "${legacy_dir}" "${framework_dir}"; do
        rm -rf "${dir}"
        mkdir -p "${dir}"

        local item base_item
        for item in "${SCRIPT_DIR}"/*; do
            base_item=$(basename "${item}")
            if [[ "${base_item}" == "regressionTests" ]]; then
                continue
            fi
            cp -a "${item}" "${dir}/"
        done

        sed -i "s|^endTime .*|endTime         ${COMPARISON_END_TIME};|" \
            "${dir}/system/controlDict"

        sed -i 's|solvePressureEqn[[:space:]]*yes;|solvePressureEqn no;|' \
            "${dir}/constant/mechanicalProperties"

        if grep -q "^writePrecision" "${dir}/system/controlDict"; then
            sed -i 's|^writePrecision.*|writePrecision  14;|' \
                "${dir}/system/controlDict"
        else
            echo "writePrecision  14;" >> "${dir}/system/controlDict"
        fi
    done

    sed -i \
        's|^\( *\)nCorrectors|\1useMechanicalConstitutiveLawManager yes;\n\1nCorrectors|' \
        "${framework_dir}/constant/solidProperties"

    for dir in "${legacy_dir}" "${framework_dir}"; do
        ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
            echo "FAIL: the comparison could not run ${dir}"
            return 1
        }
    done

    if solids4Foam::regressionCaseSkipped "${legacy_dir}/${ALLRUN_LOGFILE}"; then
        echo "Skipping the framework comparison: the case skipped here"
        return 0
    fi

    if ! grep -q "Selecting mechanical constitutive law" \
        "${framework_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the framework arm did not use the framework"
        return 1
    fi

    if grep -q "Selecting mechanical constitutive law" \
        "${legacy_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the legacy arm used the framework"
        return 1
    fi

    local t
    t=$(latest_time_dir "${legacy_dir}")

    if [[ "${t}" != "$(latest_time_dir "${framework_dir}")" ]]; then
        echo "FAIL: the arms reached different times"
        return 1
    fi

    if [[ -z "${t}" || ! -f "${legacy_dir}/${t}/D" \
       || ! -f "${framework_dir}/${t}/D" ]]
    then
        echo "FAIL: the comparison produced no D field"
        return 1
    fi

    # Round-off rather than bit-identical: the two arms reach the same answer
    # by different orderings of the same arithmetic
    local rel
    rel=$(awk '
        function abs(x) { return x < 0 ? -x : x }
        FNR == NR {
            if ($0 ~ /^\(/)
            {
                gsub(/[()]/, ""); n++
                for (i = 1; i <= NF; i++) a[n, i] = $i
            }
            next
        }
        {
            if ($0 ~ /^\(/)
            {
                gsub(/[()]/, ""); k++
                for (i = 1; i <= NF; i++)
                {
                    d = abs($i - a[k, i]); if (d > maxd) maxd = d
                    v = abs(a[k, i]);      if (v > maxv) maxv = v
                }
            }
        }
        END { printf "%.10g\n", (maxv > 0 ? maxd/maxv : maxd) }
    ' "${legacy_dir}/${t}/D" "${framework_dir}/${t}/D")

    if awk "BEGIN {exit !(${rel} < 1e-10)}"; then
        printf "PASS: framework and legacy agree to round-off, %.3g\n" "${rel}"
        return 0
    fi

    printf "FAIL: framework and legacy differ, relative D diff = %.4g\n" "${rel}"
    return 1
}

failures=0

if awk "BEGIN {exit !(${final_force_y} >= ${FORCE_Y_MIN} && ${final_force_y} <= ${FORCE_Y_MAX})}"; then
    printf "PASS: Final force_y = %.6g\n" "${final_force_y}"
else
    printf "FAIL: Final force_y = %.6g\n" "${final_force_y}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_disp_y} >= ${DISP_Y_MIN} && ${final_disp_y} <= ${DISP_Y_MAX})}"; then
    printf "PASS: Final disp_y = %.6g\n" "${final_disp_y}"
else
    printf "FAIL: Final disp_y = %.6g\n" "${final_disp_y}"
    failures=$((failures + 1))
fi

if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
fi

echo

if ! run_framework_comparison; then
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
