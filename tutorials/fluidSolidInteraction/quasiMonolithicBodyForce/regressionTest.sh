#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# quasiMonolithicBodyForce regression test
#
# hydrostatic: a uniform fluid body force g = (0 -10 0) in an enclosed box
# below a clamped solid layer. The exact discrete solution is the fluid at
# rest with the linear pressure p = g (1 - y), zero at the interface, and an
# unloaded, undeformed solid. It needs the body force in the momentum
# residual and in the wall-normal momentum balance that sets the wall and
# interface pressure gradient.
#
# gravity, bodyForce: the solid loaded by gravity (constant/solid/g) and by
# the same load as an fvOptions body force must give the same displacement.
# ============================================================

END_TIME=0.03
U_TOL=1e-9
P_TOL=1e-8
D_TOL=1e-12
D_REL_TOL=1e-8

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "quasiMonolithicBodyForce regression test"
echo "============================================================"
echo

run_variant() {
    local variant=$1
    local case_dir="${REGRESSION_ROOT}/${variant}"

    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"
    for item in "${SCRIPT_DIR}"/*; do
        if [[ "$(basename "${item}")" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done

    (cd "${case_dir}" && ./Allrun "${variant}" > "${ALLRUN_LOGFILE}" 2>&1)
}

max_abs() {
    awk '{ v = ($1 < 0 ? -$1 : $1); if (v > m) m = v } END { printf "%.6e\n", m }'
}

failed=0

for variant in hydrostatic gravity bodyForce
do
    if ! run_variant "${variant}"; then
        echo "FAIL: ${variant}: Allrun failed"
        tail -n 20 "${REGRESSION_ROOT}/${variant}/${ALLRUN_LOGFILE}" || true
        failed=1
    fi
done

if solids4Foam::regressionCaseSkipped \
    "${REGRESSION_ROOT}/hydrostatic/${ALLRUN_LOGFILE}"
then
    echo "Regression test SKIPPED"
    exit 0
fi

check() {
    local label=$1
    local value=$2
    local tol=$3

    if awk "BEGIN {exit !(${value} <= ${tol})}"; then
        echo "PASS: ${label}: ${value} <= ${tol}"
    else
        echo "FAIL: ${label}: ${value} > ${tol}"
        failed=1
    fi
}

# Every run must reach the end time and write the fields that are checked
for file in \
    hydrostatic/${END_TIME}/fluid/U hydrostatic/${END_TIME}/fluid/p \
    hydrostatic/${END_TIME}/fluid/Cy hydrostatic/${END_TIME}/solid/D \
    gravity/${END_TIME}/solid/D bodyForce/${END_TIME}/solid/D
do
    if [[ ! -f "${REGRESSION_ROOT}/${file}" ]]; then
        echo "FAIL: ${file} was not written"
        failed=1
    fi
done

if [[ ${failed} -eq 0 ]]; then
    hs="${REGRESSION_ROOT}/hydrostatic/${END_TIME}"

    check "hydrostatic: max |U| in the fluid" \
        "$(solids4Foam::internalFieldValues "${hs}/fluid/U" | max_abs)" \
        "${U_TOL}"

    # The pressure against p = 10 (1 - y) at the cell centres
    check "hydrostatic: max |p - 10 (1 - y)|" \
        "$(paste \
            <(solids4Foam::internalFieldValues "${hs}/fluid/p") \
            <(solids4Foam::internalFieldValues "${hs}/fluid/Cy") \
          | awk '{print $1 - 10*(1 - $2)}' | max_abs)" \
        "${P_TOL}"

    check "hydrostatic: max |D| in the solid" \
        "$(solids4Foam::internalFieldValues "${hs}/solid/D" | max_abs)" \
        "${D_TOL}"

    dg="${REGRESSION_ROOT}/gravity/${END_TIME}/solid/D"
    db="${REGRESSION_ROOT}/bodyForce/${END_TIME}/solid/D"
    dMax=$(solids4Foam::internalFieldValues "${dg}" | max_abs)
    dDiff=$(paste \
            <(solids4Foam::internalFieldValues "${dg}") \
            <(solids4Foam::internalFieldValues "${db}") \
          | awk '{print $1 - $2}' | max_abs)

    if awk "BEGIN {exit !(${dMax} > 0)}"; then
        check "solid: max |D(gravity) - D(bodyForce)| / max |D|" \
            "$(awk "BEGIN {printf \"%.6e\", ${dDiff}/${dMax}}")" \
            "${D_REL_TOL}"
    else
        echo "FAIL: solid: the gravity run did not deform the solid"
        failed=1
    fi
fi

echo "============================================================"
if [[ ${failed} -ne 0 ]]; then
    echo "Regression test FAILED"
    echo "============================================================"
    exit 1
fi
echo "Regression test PASSED"
echo "============================================================"
