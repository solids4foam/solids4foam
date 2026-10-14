#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
SOLIDS4FOAM_SCRIPTS="${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

if [[ -f "${SOLIDS4FOAM_SCRIPTS}" ]]; then
    source "${SOLIDS4FOAM_SCRIPTS}"
fi

# ============================================================
# heatedPlatePlaneStress regression test
# Uniform heating of a plane stress plate, free and clamped
# ============================================================

# ------------------------------------------------------------
# Analytical solution
# ------------------------------------------------------------

# E = 200 GPa, nu = 0.3, alpha = 1e-5 1/K and T - T0 = 100 K, so the free
# thermal strain is alpha*(T - T0) = 1e-3
#
# free: zero stress and D = 1e-3*(x, y, 0). On the 10 x 10 mesh of the unit
# square the cell centres are 0.05, 0.15, ..., 0.95, so the largest component
# of D is 9.5e-4 and the mean magnitude of the three components is
# (5e-4 + 5e-4 + 0)/3
#
# clamped: zero displacement and sigma_xx = sigma_yy =
# -E*alpha*(T - T0)/(1 - nu). sigma_zz is not post-processed in plane stress
# and has the same value, so the largest component of sigma is that value and
# the mean magnitude of its six components is half of it
#
# Before #137 was fixed, the thermal stress was (1 + 2*nu)/(1 + nu) = 1.23
# times too large in plane stress, and both checks failed by 23%
REF_FREE_D_MAX=9.5e-4
REF_FREE_D_MEAN=3.33333333333333e-4
REF_CLAMPED_SIGMA_MAX=2.85714285714286e8
REF_CLAMPED_SIGMA_MEAN=1.42857142857143e8
REL_TOL=1e-4

# ------------------------------------------------------------
# Log files
# ------------------------------------------------------------

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "heatedPlatePlaneStress regression test"
echo "free: D = alpha*(T - T0)*(x, y)"
echo "clamped: sigma_xx = sigma_yy = -E*alpha*(T - T0)/(1 - nu)"
echo "============================================================"
echo

prepare_case() {
    local case_dir=$1

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done
}

run_case() {
    local case_dir=$1

    ( cd "${case_dir}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${case_dir}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
}

# Print the latest time of a case, or fail if it did not reach the end time
latest_time() {
    local case_dir=$1
    local t end_time

    if ! grep -q "Selecting mechanical constitutive law" \
        "${case_dir}/${SOLVER_LOGFILE}"
    then
        echo "FAIL: ${case_dir##*/}: no mechanical constitutive law" >&2
        return 1
    fi

    t=$(solids4Foam::latestTime "${case_dir}")
    end_time=$(sed -n 's/^endTime[[:space:]]*\([^;]*\);.*/\1/p' \
        "${case_dir}/system/controlDict")

    if [[ -z "${t}" || -z "${end_time}" ]] \
        || ! awk "BEGIN {exit !((${t} - ${end_time})^2 <= 1e-20)}"
    then
        echo "FAIL: ${case_dir##*/}: stopped at '${t}', not at '${end_time}'" >&2
        return 1
    fi

    echo "${t}"
}

# The norms are of magnitudes, so the sign is checked separately: the signed
# mean of the components is the mean magnitude, with the sign of the solution
check_signed_mean() {
    local label=$1
    local field_file=$2
    local ref_mean=$3
    local values mean

    if ! values=$(solids4Foam::internalFieldValues "${field_file}"); then
        echo "FAIL: ${label}: no field to check the sign of"
        return 1
    fi

    mean=$(awk '{ sum += $1; n++ } END { printf "%.15g\n", sum/n }' \
        <<< "${values}")

    if awk "BEGIN {
            d = ${mean} - (${ref_mean}); if (d < 0) d = -d
            r = ${ref_mean}; if (r < 0) r = -r
            exit !(d <= ${REL_TOL}*r)
        }"
    then
        printf "PASS: %s signed mean %.15g (%.15g)\n" \
            "${label}" "${mean}" "${ref_mean}"
        return 0
    fi

    printf "FAIL: %s signed mean %.15g (%.15g)\n" \
        "${label}" "${mean}" "${ref_mean}"
    return 1
}

failures=0

# ------------------------------------------------------------
# Free plate: the tutorial as it is
# ------------------------------------------------------------

FREE_DIR="${REGRESSION_ROOT}/free"
prepare_case "${FREE_DIR}"
run_case "${FREE_DIR}"

if t=$(latest_time "${FREE_DIR}"); then
    solids4Foam::checkFieldNorms "free: D" "${FREE_DIR}/${t}/D" \
        "${REF_FREE_D_MAX}" "${REF_FREE_D_MEAN}" "${REL_TOL}" \
        || failures=$((failures + 1))
    check_signed_mean "free: D" "${FREE_DIR}/${t}/D" "${REF_FREE_D_MEAN}" \
        || failures=$((failures + 1))
else
    failures=$((failures + 1))
fi

( cd "${FREE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

# ------------------------------------------------------------
# Clamped plate: the traction-free edges are fixed instead
# ------------------------------------------------------------

CLAMPED_DIR="${REGRESSION_ROOT}/clamped"
prepare_case "${CLAMPED_DIR}"
sed 's/solidTraction;/fixedDisplacement;/' "${CLAMPED_DIR}/0/D" \
    > "${CLAMPED_DIR}/0/D.clamped"
mv "${CLAMPED_DIR}/0/D.clamped" "${CLAMPED_DIR}/0/D"
run_case "${CLAMPED_DIR}"

if t=$(latest_time "${CLAMPED_DIR}"); then
    solids4Foam::checkFieldNorms "clamped: sigma" "${CLAMPED_DIR}/${t}/sigma" \
        "${REF_CLAMPED_SIGMA_MAX}" "${REF_CLAMPED_SIGMA_MEAN}" "${REL_TOL}" \
        || failures=$((failures + 1))
    check_signed_mean "clamped: sigma" "${CLAMPED_DIR}/${t}/sigma" \
        "-${REF_CLAMPED_SIGMA_MEAN}" || failures=$((failures + 1))
else
    failures=$((failures + 1))
fi

( cd "${CLAMPED_DIR}" && ./Allclean > /dev/null 2>&1 ) || true

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

echo
if (( failures == 0 ))
then
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
