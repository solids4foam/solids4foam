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
# multiMaterialPlate regression test
# Runs the default case, with a mid-layer of thickness
# h = 0.01W, and checks the normal stress sampled along the
# two interfaces.
# ============================================================

W=0.03
SIGMA0=1e6
E1=206e9
E2=70.3e9
E3=4.93e9
NU1=0.3
NU2=0.345
NU3=0.33

# Equilibrium: the mean of the normal stress over an interface is the applied
# stress, whatever the materials. It is 1 +- 2e-6 on OpenFOAM.com v2412
MEAN_STRESS_ERR_MAX=1e-3

# The order of the singularity is fitted in R_MIN < r/W < R_MAX and compared
# with the analytical order. Closer to the free edge than R_MIN, the stress is
# not resolved by the mesh. The singular term of interface 23 is not yet as
# dominant in this range as the one of interface 13: on OpenFOAM.com v2412 the
# fitted orders are 0.2876 and 0.2191, against 0.2823 and 0.2392
R_MIN=1e-4
R_MAX=3e-4
ORDER13_REL_ERR_MAX=0.05
ORDER23_REL_ERR_MAX=0.12

# Further from the free edge, the stress has a local minimum and then a local
# maximum, where the two interfaces carry nearly the same stress. These are
# the values of OpenFOAM.com v2412, which agree with Figure 11 of Tuković et
# al. (2013)
REF_MIN13=0.99869
REF_MAX13=1.09742
REF_MIN23=0.80732
REF_MAX23=1.09989
EXTREMA_REL_ERR_MAX=0.01

SOLVER_LOGFILE="log.solids4Foam"
ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "multiMaterialPlate regression test"
echo "Mean interface stress error   < ${MEAN_STRESS_ERR_MAX}"
echo "Singularity order 13 error    < ${ORDER13_REL_ERR_MAX}"
echo "Singularity order 23 error    < ${ORDER23_REL_ERR_MAX}"
echo "Local stress extrema error    < ${EXTREMA_REL_ERR_MAX}"
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
    find "${CASE_DIR}" -name "$1_sigma.xy" | sort | tail -n 1
}

# The order p of the stress singularity, sigma ~ r^(-p), where the interface
# between two materials meets a free edge, in plane strain: the root
# lambda = 1 - p of the characteristic equation of Bogy (1970), found by
# bisection between lambda = 0.5 and the trivial root lambda = 1
analytical_order() {
    awk -v Ea="$1" -v nua="$2" -v Eb="$3" -v nub="$4" '
    function f(l,    s) {
        s = sin(pi*l/2)^2
        return l*l*(l*l - 1)*a*a + 2*l*l*(s - l*l)*a*b \
            + (s - l*l)^2*b*b + s*(1 - s)
    }
    BEGIN {
        pi = atan2(0, -1)
        k = (Eb/(2*(1 + nub)))/(Ea/(2*(1 + nua)))
        ma = 4*(1 - nua)
        mb = 4*(1 - nub)
        a = (mb - k*ma)/(mb + k*ma)
        b = ((mb - 2) - k*(ma - 2))/(mb + k*ma)
        lo = 0.5
        hi = 0.9999
        for (i = 0; i < 60; i++) {
            mid = 0.5*(lo + hi)
            if (f(lo)*f(mid) <= 0) { hi = mid } else { lo = mid }
        }
        printf "%.6f\n", 1 - 0.5*(lo + hi)
    }'
}

# From a sampled interface stress, as "mean<TAB>order<TAB>min<TAB>max": the
# mean of sigma_xx/sigma0 over the interface, the order fitted in
# R_MIN < r/W < R_MAX, the minimum in 1e-4 < r/W < 1e-2 and the maximum in
# 2e-3 < r/W < 1e-1
interface_stress() {
    awk -v W="${W}" -v sigma0="${SIGMA0}" -v rMin="${R_MIN}" -v rMax="${R_MAX}" '
    NF >= 7 {
        r = ($1 + W)/W
        s = $2/sigma0

        if (n > 0) { integral += 0.5*(s + sPrev)*(r - rPrev) }
        n++
        rPrev = r
        sPrev = s

        if (r > rMin && r < rMax && s > 0) {
            x = log(r); y = log(s)
            m++; sx += x; sy += y; sxx += x*x; sxy += x*y
        }
        if (r > 1e-4 && r < 1e-2 && (!haveMin || s < sMin)) {
            sMin = s; haveMin = 1
        }
        if (r > 2e-3 && r < 1e-1 && (!haveMax || s > sMax)) {
            sMax = s; haveMax = 1
        }
    }
    END {
        if (n < 2 || m < 2 || !haveMin || !haveMax) { exit 1 }
        printf "%.8g\t%.8g\t%.8g\t%.8g\n", integral, \
            -(m*sxy - sx*sy)/(m*sxx - sx*sx), sMin, sMax
    }' "$1"
}

# PASS if |value - reference| < tolerance*|scale|
check_value() {
    local label="$1" value="$2" reference="$3" scale="$4" tolerance="$5"
    local err

    err=$(awk "BEGIN {d = (${value}) - (${reference}); \
        if (d < 0) d = -d; s = (${scale}); if (s < 0) s = -s; \
        printf \"%.6g\", d/s}")

    if awk "BEGIN {exit !(${err} < ${tolerance})}"; then
        printf "PASS: %s = %.6g (reference %.6g, error %.3g)\n" \
            "${label}" "${value}" "${reference}" "${err}"
    else
        printf "FAIL: %s = %.6g (reference %.6g, error %.3g)\n" \
            "${label}" "${value}" "${reference}" "${err}"
        failures=$((failures + 1))
    fi
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

failures=0

# Three materials, each with its own law
n_laws=$({ grep -c "Selecting mechanical constitutive law" \
    "${CASE_DIR}/${SOLVER_LOGFILE}" 2> /dev/null || true; })

if [[ "${n_laws}" == 3 ]]; then
    echo "PASS: the case constructed three mechanical constitutive laws"
else
    echo "FAIL: the case constructed ${n_laws:-0} mechanical constitutive" \
        "laws, not three"
    failures=$((failures + 1))
fi

for interface in 13 23; do
    data_file="$(sample_file "interface${interface}")"

    if [[ -z "${data_file}" || ! -f "${data_file}" ]]; then
        echo "FAIL: Could not find the sampled stress of interface ${interface}"
        failures=$((failures + 1))
        continue
    fi

    if ! result=$(interface_stress "${data_file}"); then
        echo "FAIL: Could not extract the stress of interface ${interface}"
        failures=$((failures + 1))
        continue
    fi

    IFS=$'\t' read -r mean order s_min s_max <<< "${result}"

    if [[ "${interface}" == 13 ]]; then
        p=$(analytical_order "${E1}" "${NU1}" "${E3}" "${NU3}")
        order_tol="${ORDER13_REL_ERR_MAX}"
        ref_min="${REF_MIN13}"
        ref_max="${REF_MAX13}"
    else
        p=$(analytical_order "${E2}" "${NU2}" "${E3}" "${NU3}")
        order_tol="${ORDER23_REL_ERR_MAX}"
        ref_min="${REF_MIN23}"
        ref_max="${REF_MAX23}"
    fi

    check_value "Interface ${interface}: mean stress" \
        "${mean}" 1 1 "${MEAN_STRESS_ERR_MAX}"
    check_value "Interface ${interface}: singularity order" \
        "${order}" "${p}" "${p}" "${order_tol}"
    check_value "Interface ${interface}: local minimum stress" \
        "${s_min}" "${ref_min}" "${ref_min}" "${EXTREMA_REL_ERR_MAX}"
    check_value "Interface ${interface}: local maximum stress" \
        "${s_max}" "${ref_max}" "${ref_max}" "${EXTREMA_REL_ERR_MAX}"
done

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
