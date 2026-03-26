#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# slidingFrictionBall regression test
# Uses a short point-displacement history as a cheap contact check.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

REG_END_TIME=5
DX_MIN=0.00199
DX_MAX=0.00201

ALLRUN_LOGFILE="log.Allrun"
POINT_FILE="postProcessing/0/solidPointDisplacement_pointDisp.dat"

echo "============================================================"
echo "slidingFrictionBall regression test"
echo "Regression end time = ${REG_END_TIME}"
echo "Final Dx in [${DX_MIN}, ${DX_MAX}] m"
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

shorten_case() {
    local controlDict="${CASE_DIR}/system/controlDict"
    sed -i.bak 's/^endTime[[:space:]]\+100;/endTime         5;/' "${controlDict}"
    rm -f "${controlDict}.bak"

    cat >> "${controlDict}" <<'EOF'

functions
{
    pointDisp
    {
        type    solidPointDisplacement;
        point   (-1.2 0 0);
    }
}
EOF
}

# ------------------------------------------------------------
# Clean & run case
# ------------------------------------------------------------

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
    shorten_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

if [[ ! -f "${CASE_DIR}/${POINT_FILE}" ]]; then
    echo "FAIL: Could not find ${POINT_FILE}"
    exit 1
fi

final_dx=$(awk 'END {print $2}' "${CASE_DIR}/${POINT_FILE}")

if [[ -z "${final_dx}" ]]; then
    echo "FAIL: Could not extract final Dx"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${final_dx} >= ${DX_MIN} && ${final_dx} <= ${DX_MAX})}"; then
    printf "PASS: Final Dx = %.6g\n" "${final_dx}"
else
    printf "FAIL: Final Dx = %.6g\n" "${final_dx}"
    failures=$((failures + 1))
fi

# Clean case again
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
