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
# layeredPipe regression test
# Compares sampled transformed stresses against the analytical
# cylinder solution used by the tutorial plots.
# ============================================================

RADIUS_STRESS_ERR_MAX=0.03
THETA_POINT_ERR_MAX=0.01

R1=0.05
R2=0.07
R3=0.1
E1=20e9
E2=200e9
NU1=0.35
NU2=0.3
PLOAD=1e5

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"
CONSTITUTIVE_LOGFILE="log.Test-mechanicalConstitutiveLaw"

echo "============================================================"
echo "layeredPipe regression test"
echo "Max radial stress error  < ${RADIUS_STRESS_ERR_MAX}"
echo "Outer-point theta error  < ${THETA_POINT_ERR_MAX}"
echo "Plus the mechanicalConstitutiveLaw framework checks"
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
}

sample_file() {
    local preferred
    preferred="$(find "${CASE_DIR}" -name 'line_sigma:Transformed.xy' | sort | tail -n 1)"
    if [[ -n "${preferred}" ]]; then
        printf '%s\n' "${preferred}"
        return 0
    fi

    find "${CASE_DIR}" -name 'sigma:Transformed' | sort | tail -n 1
}

# Exercise the mechanicalConstitutiveLawManager evaluation paths on this case.
# This tutorial is used because it is the only one with more than one material,
# and the manager's integration-point addressing is per material. The framework
# is not yet used by any solid model, so this is its only runtime coverage
run_constitutive_test() {
    if ! command -v Test-mechanicalConstitutiveLaw > /dev/null 2>&1; then
        echo "SKIP: Test-mechanicalConstitutiveLaw not found in PATH"
        return 0
    fi

    if [[ ! -d "${CASE_DIR}/constant/polyMesh" ]]; then
        echo "SKIP: mechanicalConstitutiveLaw checks (case has no mesh)"
        return 0
    fi

    if ( cd "${CASE_DIR}" && Test-mechanicalConstitutiveLaw \
            > "${CONSTITUTIVE_LOGFILE}" 2>&1 )
    then
        local n_passed
        n_passed=$(grep -c 'PASS:' "${CASE_DIR}/${CONSTITUTIVE_LOGFILE}" || true)
        echo "PASS: mechanicalConstitutiveLaw checks (${n_passed} checks)"
        return 0
    fi

    echo "FAIL: mechanicalConstitutiveLaw checks"
    grep 'FAIL:' "${CASE_DIR}/${CONSTITUTIVE_LOGFILE}" || true
    return 1
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

DATA_FILE="$(sample_file)"

if [[ -z "${DATA_FILE}" || ! -f "${DATA_FILE}" ]]; then
    echo "FAIL: Could not find sampled stress data"
    exit 1
fi

max_radial_err="$(awk -v r1="${R1}" -v r2="${R2}" -v r3="${R3}" \
        -v e1="${E1}" -v e2="${E2}" -v nu1="${NU1}" -v nu2="${NU2}" -v p="${PLOAD}" '
    function abs(x) { return x < 0 ? -x : x }
    BEGIN {
        pint = (2*r1*r1*p/(e1*(r2*r2-r1*r1))) / (((1.0/e2)*(((r3*r3+r2*r2)/(r3*r3-r2*r2))+nu2)) + ((1.0/e1)*(((r2*r2+r1*r1)/(r2*r2-r1*r1))-nu1)))
        maxRad = 0
    }
    NF >= 5 {
        r = $1
        sigmaR = $2

        if (r < r2) {
            sigmaRAnal = (r1*r1*p - r2*r2*pint + (pint-p)*(r1*r2/r)^2) / (r2*r2-r1*r1)
        } else {
            sigmaRAnal = (r2*r2*pint - pint*(r2*r3/r)^2) / (r3*r3-r2*r2)
        }

        radialErr = abs(sigmaR - sigmaRAnal) / p
        if (radialErr > maxRad) {
            maxRad = radialErr
        }
    }
    END {
        printf "%.12g\n", maxRad
    }
    ' "${DATA_FILE}")"

outer_theta_err="$(awk -v r1="${R1}" -v r2="${R2}" -v r3="${R3}" \
        -v e1="${E1}" -v e2="${E2}" -v nu1="${NU1}" -v nu2="${NU2}" -v p="${PLOAD}" '
    function abs(x) { return x < 0 ? -x : x }
    BEGIN {
        pint = (2*r1*r1*p/(e1*(r2*r2-r1*r1))) / (((1.0/e2)*(((r3*r3+r2*r2)/(r3*r3-r2*r2))+nu2)) + ((1.0/e1)*(((r2*r2+r1*r1)/(r2*r2-r1*r1))-nu1)))
        lastThetaAnal = 0
        lastThetaSample = 0
    }
    NF >= 5 {
        r = $1
        sigmaTheta = $5

        if (r < r2) {
            sigmaThetaAnal = (r1*r1*p - r2*r2*pint - (pint-p)*(r1*r2/r)^2) / (r2*r2-r1*r1)
        } else {
            sigmaThetaAnal = (r2*r2*pint + pint*(r2*r3/r)^2) / (r3*r3-r2*r2)
        }

        lastThetaAnal = sigmaThetaAnal
        lastThetaSample = sigmaTheta
    }
    END {
        thetaErr = abs(lastThetaSample - lastThetaAnal) / p
        printf "%.12g\n", thetaErr
    }
    ' "${DATA_FILE}")"

if [[ -z "${max_radial_err}" || -z "${outer_theta_err}" ]]; then
    echo "FAIL: Could not extract regression errors"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${max_radial_err} < ${RADIUS_STRESS_ERR_MAX})}"; then
    printf "PASS: Max radial stress error = %.6g\n" "${max_radial_err}"
else
    printf "FAIL: Max radial stress error = %.6g\n" "${max_radial_err}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${outer_theta_err} < ${THETA_POINT_ERR_MAX})}"; then
    printf "PASS: Outer-point theta error = %.6g\n" "${outer_theta_err}"
else
    printf "FAIL: Outer-point theta error = %.6g\n" "${outer_theta_err}"
    failures=$((failures + 1))
fi

if ! run_constitutive_test; then
    failures=$((failures + 1))
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
