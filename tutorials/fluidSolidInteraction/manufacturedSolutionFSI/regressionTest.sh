#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# ============================================================
# manufacturedSolutionFSI regression test
# ============================================================

# The coarse case (16 cells along the interface, deltaT = 0.02) is run with the
# Aitken and Robin-Neumann couplings, and the L2 errors against the
# manufactured solution at t = 0.4 are compared with the development
# reference errors of
# verification/reference/manufacturedSolutionFSI_verification_references.json.
# The solid uses tractionPatchSurfaceGradient, which is ignored until
# solids4foam PR #756 is merged; these references are then to be updated
REL_TOL=0.02

COUPLINGS=(aitken robin)

# L2 errors at t = 0.4: U p D intTrac
ref_errors() {
    case "$1" in
        aitken) echo "0.038365 0.06739 0.0037206 0.088893" ;;
        robin) echo "0.043958 0.086565 0.0038143 0.11443" ;;
    esac
}

ALLRUN_LOGFILE="log.Allrun"
ERRORS_FILE="postProcessing/mmsFsiErrors/errors.dat"

CHECK_ONLY=false
for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run) CHECK_ONLY=true ;;
    esac
done

prepare_case() {
    local case_dir="$1"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"
    for item in "${SCRIPT_DIR}"/*; do
        name=$(basename "${item}")
        case "${name}" in
            regressionTests|verification|derivation|src) continue ;;
            postProcessing|log.*|[1-9]*) continue ;;
        esac
        cp -a "${item}" "${case_dir}/"
    done
    rm -rf "${case_dir}"/constant/*/polyMesh

    # Share the tutorial library, so that it is compiled only once
    ln -s "${SCRIPT_DIR}/src" "${case_dir}/src"
}

failures=0
for coupling in "${COUPLINGS[@]}"; do
    case_dir="${REGRESSION_ROOT}/${coupling}"

    if [[ "${CHECK_ONLY}" == false ]]; then
        prepare_case "${case_dir}"
        (cd "${case_dir}" && ./Allrun "${coupling}" > "${ALLRUN_LOGFILE}" 2>&1)
    fi

    # A skip is only valid if the tutorial declared one in the Allrun log
    if solids4Foam::regressionCaseSkipped "${case_dir}/${ALLRUN_LOGFILE}"; then
        echo "Skipping regression checks because the tutorial skipped in this environment"
        exit 0
    fi

    solver_log="${case_dir}/log.solids4Foam"
    if [[ ! -f "${solver_log}" ]] || ! grep -q '^End' "${solver_log}"; then
        echo "FAIL: missing or incomplete solver log: ${solver_log}"
        failures=$((failures + 1))
        continue
    fi

    # Columns of errors.dat: time, then U_L2 is 3, p_L2 6, D_L2 10 and
    # intTrac_L2 19
    errors=$(awk '$1 == 0.4 {print $3, $6, $10, $19}' \
        "${case_dir}/${ERRORS_FILE}")
    if [[ -z "${errors}" ]]; then
        echo "FAIL: no errors at t = 0.4 in ${case_dir}/${ERRORS_FILE}"
        failures=$((failures + 1))
        continue
    fi

    echo
    echo "Checking ${coupling}"
    IFS=' ' read -r -a values <<< "${errors}"
    IFS=' ' read -r -a refs <<< "$(ref_errors "${coupling}")"
    names=(U_L2 p_L2 D_L2 intTrac_L2)
    for i in "${!names[@]}"; do
        if awk -v v="${values[$i]}" -v r="${refs[$i]}" -v tol="${REL_TOL}" \
            'BEGIN {d = v - r; if (d < 0) d = -d; exit !(v > 0 && d <= tol*r)}'
        then
            printf 'PASS: %s = %.6g (reference %.6g)\n' \
                "${names[$i]}" "${values[$i]}" "${refs[$i]}"
        else
            printf 'FAIL: %s = %.6g (reference %.6g, relative tolerance %g)\n' \
                "${names[$i]}" "${values[$i]}" "${refs[$i]}" "${REL_TOL}"
            failures=$((failures + 1))
        fi
    done
done

echo
if (( failures == 0 )); then
    echo "Regression test PASSED"
    exit 0
else
    echo "Regression test FAILED (${failures} checks)"
    exit 1
fi
