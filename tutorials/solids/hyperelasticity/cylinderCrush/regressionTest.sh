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

# Against the removed legacy mechanicalModel.
#
# The comparison run turns the pressure equation off. The case ships with
# solvePressureEqn yes and a smoothing scale factor; the legacy law smoothed the
# hydrostatic stress itself, and on the framework the solid model does it
# instead. Turning it off here checks the MooneyRivlin law on its own. The
# shipped configuration, smoothing included, is covered by the force band
# above - without the smoothing the final force is about -497 N, outside it -
# and by the log line checked below.
#
# Two of the thirty steps. The legacy solver reached its corrector limit from
# the third step onwards, so past that it was an unconverged answer
COMPARISON_END_TIME=2

# The norms of the comparison run's final D, max and mean component magnitude,
# on the removed legacy mechanicalModel, from the last commit that had it
# (mcl-stage8-coverage, c3a92b3d), on foam-extend 4.1, the one fork this case
# runs on. The framework reproduced the legacy D field to round-off, 1e-10 of
# its largest value, and that is the tolerance here
LEGACY_D_MAX=0.0066666666744038
LEGACY_D_MEAN=0.00116545530838457
LEGACY_D_REL_TOL=1e-10

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

check_completed() {
    local case_dir="$1"
    local expected_time="$2"
    local actual_time

    if ! grep -q "^End" "${case_dir}/${SOLVER_LOGFILE}" \
      || grep -qE "Nonlinear solve did not converge|SNES convergence error|FOAM FATAL" \
          "${case_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${case_dir} did not complete and converge"
        return 1
    fi

    actual_time=$(latest_time_dir "${case_dir}")
    if ! awk "BEGIN {exit !((${actual_time:-0} - ${expected_time})^2 <= 1e-20)}"
    then
        echo "FAIL: ${case_dir} stopped at ${actual_time:-none}; expected ${expected_time}"
        return 1
    fi
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

run_legacy_comparison() {
    local dir="${REGRESSION_ROOT}/comparison"

    if [ "$CHECK_ONLY" = false ]; then
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

        ( cd "${dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || {
            echo "FAIL: the comparison could not run ${dir}"
            return 1
        }
    else
        echo "Check-only: checking the existing comparison run"
    fi

    if solids4Foam::regressionCaseSkipped "${dir}/${ALLRUN_LOGFILE}"; then
        echo "Skipping the legacy comparison: the case skipped here"
        return 0
    fi

    if ! grep -q "Selecting mechanical constitutive law" \
        "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the comparison run constructed no mechanical constitutive law"
        return 1
    fi

    # Without solvePressureEqn there must be no smoothing, or this is not the
    # law on its own
    if grep -q "smoothing the hydrostatic stress" "${dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: the comparison run smoothed the hydrostatic stress"
        return 1
    fi

    check_completed "${dir}" "${COMPARISON_END_TIME}" || return 1

    check_field_against_legacy "D at t = ${COMPARISON_END_TIME}" \
        "${dir}/$(latest_time_dir "${dir}")/D" \
        "${LEGACY_D_MAX}" "${LEGACY_D_MEAN}" "${LEGACY_D_REL_TOL}"
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

# The shipped case asks for the hydrostatic stress smoothing, and the solid
# model has to say it is doing it
if grep -q "smoothing the hydrostatic stress (solvePressureEqn)" \
    "${CASE_DIR}/${SOLVER_LOGFILE}"
then
    echo "PASS: the hydrostatic stress is smoothed"
else
    echo "FAIL: no hydrostatic stress smoothing in the log"
    failures=$((failures + 1))
fi

if [ "$CHECK_ONLY" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
fi

echo

if ! run_legacy_comparison; then
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
