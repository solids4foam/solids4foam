#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../../applications/scripts/solids4FoamScripts.sh"

# ============================================================
# leftVentricle regression test
#
# The full case is three cardiac cycles (30,000 steps), so the test runs the
# first 0.01 s (100 steps, the wall moving from about t = 0.003 s) on the
# coarsest mesh of the study (M1, 39,433 cells) on 4 ranks, for
# 1. the implicit coupling of order 2 (the tutorial's set-up), and
# 2. the explicit coupling of order 2,
# and checks for each:
#   - that volume is conserved between the 3D and 0D models: the change in
#     the mesh volume against the change in the circulation's volume;
#   - the mitral and aortic patch pressures at t = 0.01 s against reference
#     values.
# ============================================================

# ------------------------------------------------------------
# Regression tolerances and reference values
# ------------------------------------------------------------

END_TIME=0.01
NPROCS=4

# Mesh M1 (system/meshDict)
MAX_CELL_SIZE=1.96e-3
FIRST_LAYER=2.94e-4

# Largest |dV_LV + dV_0D| over the run (mL): 2.1e-6 (implicit) and 6.4e-6
# (explicit) against a volume change of 0.0135 mL by t = 0.01 s; the full case
# reaches 3e-5 mL over three cycles
VOLUME_TOL=1e-5

# Patch pressures at t = END_TIME (Pa) and their tolerance. The scatter with
# another platform is about 1e-3 Pa (OpenFOAM-v2412 on macOS/clang against
# Linux/GCC, 32 ranks); the two couplings differ by about 0.1 Pa
REF_I2_PMV=4977.79361832
REF_I2_PAV=5012.03396212
REF_E2_PMV=4977.88299068
REF_E2_PAV=5012.16045415
PRESSURE_TOL=0.02

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
COUPLING_FILE="postProcessing/circulationCoupling/0/coupling.dat"
VOLUME_FILE="volumeConservation.dat"

echo "============================================================"
echo "leftVentricle regression test (M1, t <= ${END_TIME} s, ${NPROCS} ranks)"
echo "|dV_LV + dV_0D| < ${VOLUME_TOL} mL, patch pressures within ${PRESSURE_TOL} Pa"
echo "============================================================"
echo

prepare_case() {
    local d="$1"
    rm -rf "${d}"
    mkdir -p "${d}"

    # Copy the inputs only: no results of a previous run of the tutorial
    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        case "${base_item}" in
            regressionTests|postProcessing|processor*|log.*|images|*.png|*.dat)
                continue
                ;;
        esac
        if [[ "${base_item}" != "0" \
           && "${base_item}" =~ ^[0-9.eE+-]+$ ]]; then
            continue
        fi
        cp -a "${item}" "${d}/"
    done
    rm -rf "${d}/constant/polyMesh"

    set_entry "${d}/system/controlDict" endTime "${END_TIME}"
    set_entry "${d}/system/controlDict" writeInterval "${END_TIME}"
    set_entry "${d}/system/decomposeParDict" numberOfSubdomains "${NPROCS}"
    set_entry "${d}/system/meshDict" maxCellSize "${MAX_CELL_SIZE}"
    set_entry "${d}/system/meshDict" boundaryLayers/maxFirstLayerThickness \
        "${FIRST_LAYER}"
    for patch in INLET OUTLET; do
        set_entry "${d}/system/meshDict" \
            "boundaryLayers/patchBoundaryLayers/${patch}/maxFirstLayerThickness" \
            "${FIRST_LAYER}"
    done
}

set_entry() {
    foamDictionary -entry "$2" -set "$3" "$1" > /dev/null
}

# ------------------------------------------------------------
# Clean & run cases
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

CASES=(implicit explicit)

if [ "$CHECK_ONLY" = false ]; then
    for scheme in "${CASES[@]}"; do
        d="${REGRESSION_ROOT}/${scheme}"
        prepare_case "${d}"
        set_entry "${d}/constant/circulationProperties" coupling/scheme \
            "${scheme}"
        ( cd "${d}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 ) || true
    done
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if solids4Foam::regressionCaseSkipped \
    "${REGRESSION_ROOT}/implicit/${ALLRUN_LOGFILE}"
then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

# The last row of a coupling.dat at t = END_TIME: "<Pmv> <Pav>", or nothing
end_pressures() {
    [[ -f "$1" ]] || return 0
    awk -v tend="${END_TIME}" '
        /^# time/ {
            for (i = 3; i <= NF; i++) {
                if ($i == "LV.inlet.Ppatch") cmv = i - 1
                if ($i == "LV.outlet.Ppatch") cav = i - 1
            }
        }
        /^#/ || NF == 0 { next }
        { d = $1 - tend; if (d < 0) d = -d; if (d < 1e-9) { pmv = $cmv; pav = $cav } }
        END { if (cmv > 0 && cav > 0 && pmv != "") print pmv, pav }
    ' "$1"
}

# check_pressure <label> <value> <reference>
check_pressure() {
    if awk -v v="$2" -v r="$3" -v tol="${PRESSURE_TOL}" \
        'BEGIN { d = v - r; exit !((d < 0 ? -d : d) < tol) }'
    then
        printf "PASS: %s = %s Pa (reference %s)\n" "$1" "$2" "$3"
    else
        printf "FAIL: %s = %s Pa (reference %s)\n" "$1" "$2" "$3"
        failures=$((failures + 1))
    fi
}

for scheme in "${CASES[@]}"; do
    d="${REGRESSION_ROOT}/${scheme}"

    if [[ "${scheme}" == implicit ]]; then
        refPmv=${REF_I2_PMV}
        refPav=${REF_I2_PAV}
    else
        refPmv=${REF_E2_PMV}
        refPav=${REF_E2_PAV}
    fi

    # Volume conservation over the run (written by Allrun)
    if [[ -s "${d}/${VOLUME_FILE}" ]]; then
        err=$(awk '!/^#/ { s = $4 < 0 ? -$4 : $4; if (s > m) m = s; n++ }
                   END { if (n > 0) printf "%.3g", m }' "${d}/${VOLUME_FILE}")
        if [[ -n "${err}" ]] && awk -v e="${err}" -v tol="${VOLUME_TOL}" \
            'BEGIN { exit !(e < tol) }'
        then
            printf "PASS: %s: max |dV_LV + dV_0D| = %s mL\n" "${scheme}" "${err}"
        else
            printf "FAIL: %s: max |dV_LV + dV_0D| = %s mL\n" "${scheme}" "${err}"
            failures=$((failures + 1))
        fi
    else
        echo "FAIL: ${scheme}: no ${VOLUME_FILE} (see ${d}/${ALLRUN_LOGFILE})"
        failures=$((failures + 1))
    fi

    # The patch pressures at the end of the run
    p=$(end_pressures "${d}/${COUPLING_FILE}")
    if [[ -z "${p}" ]]; then
        echo "FAIL: ${scheme}: no row at t = ${END_TIME} in ${COUPLING_FILE}"
        failures=$((failures + 1))
    else
        check_pressure "${scheme}: mitral pressure" "${p%% *}" "${refPmv}"
        check_pressure "${scheme}: aortic pressure" "${p##* }" "${refPav}"
    fi
done

# Clean cases again, keeping them when a check failed
if [[ "$CHECK_ONLY" = false && ${failures} -eq 0 ]]; then
    for scheme in "${CASES[@]}"; do
        ( cd "${REGRESSION_ROOT}/${scheme}" && ./Allclean > /dev/null 2>&1 ) \
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
