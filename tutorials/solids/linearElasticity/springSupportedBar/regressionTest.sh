#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# ============================================================
# springSupportedBar regression test
# Compares the spring-end and loaded-end displacements with the exact
# solutions, for the linearGeometry and totalLagrangian cases
# ============================================================

# Relative tolerance against the exact solutions
REL_TOL=1e-6

# Exact solutions (see README.md)
# linearGeometry: t = 1 kPa, k = 1e7 Pa/m, L = 10 mm, E = 100 kPa, nu = 0.3
LG_SPRING_END=1e-4
LG_LOADED_END=1.742857142857e-4
# totalLagrangian: F = 9 x 0.02 N, kN = 1e8 Pa/m, A0 = 9 mm^2
TL_SPRING_END=2e-4

ALLRUN_LOGFILE="log.Allrun"

echo "============================================================"
echo "springSupportedBar regression test"
echo "Relative tolerance: ${REL_TOL}"
echo "============================================================"
echo

prepare_case() {
    local case_dir="$1"
    rm -rf "${case_dir}"
    mkdir -p "${case_dir}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${case_dir}/"
    done
    # never start a variant from another variant's mesh
    \rm -rf "${case_dir}/constant/polyMesh"
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

for option in linearGeometry totalLagrangian; do
    if [ "$CHECK_ONLY" = false ]; then
        prepare_case "${REGRESSION_ROOT}/${option}"
        ( cd "${REGRESSION_ROOT}/${option}" \
            && ./Allrun "${option}" > "${ALLRUN_LOGFILE}" 2>&1 ) || true
    fi
done

# Average x displacement of a patch at the last time
# (column 8 of the solidDisplacements file)
extract_avg_dx() {
    local file="${REGRESSION_ROOT}/$1/postProcessing/0/solidDisplacements$2.dat"
    if [[ -f "${file}" ]]; then
        grep -v "^#" "${file}" | tail -n 1 | awk '{print $8}'
    fi
}

failures=0

check() {
    local name="$1" value="$2" exact="$3"

    if [[ -z "${value}" ]]; then
        echo "FAIL: ${name}: could not extract the value"
        failures=$((failures + 1))
        return
    fi

    if awk "BEGIN {e = (${value} - ${exact})/${exact}; if (e < 0) e = -e; \
        exit !(e <= ${REL_TOL})}"; then
        printf "PASS: %s = %.10g (exact %.10g)\n" "${name}" "${value}" "${exact}"
    else
        printf "FAIL: %s = %.10g (exact %.10g)\n" "${name}" "${value}" "${exact}"
        failures=$((failures + 1))
    fi
}

check "linearGeometry Dx(springEnd)" \
    "$(extract_avg_dx linearGeometry springEnd)" "${LG_SPRING_END}"
check "linearGeometry Dx(loadedEnd)" \
    "$(extract_avg_dx linearGeometry loadedEnd)" "${LG_LOADED_END}"
check "totalLagrangian Dx(springEnd)" \
    "$(extract_avg_dx totalLagrangian springEnd)" "${TL_SPRING_END}"

# Clean cases again
if [ "$CHECK_ONLY" = false ]; then
    for option in linearGeometry totalLagrangian; do
        ( cd "${REGRESSION_ROOT}/${option}" && ./Allclean > /dev/null 2>&1 ) \
            || true
    done
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
