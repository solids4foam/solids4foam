#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# Source solids4Foam scripts
source "${SCRIPT_DIR}/../../../applications/scripts/solids4FoamScripts.sh"

# GNU sed, for the in-place edits below
solids4Foam::requireGnuSed

# ============================================================
# pulmonaryArteries prestressed FSI regression test
#
# The full tutorial runs 15 prestress passes and five cardiac cycles, which
# takes about a day, so the regression test runs a copy of the case with:
#   - the tutorial meshes, unchanged: coarser meshes degrade on this geometry
#   - 3 prestress passes: the ramped first pass and two single-step passes
#   - the first few time steps of the fluid-solid interaction run
#   - 4 processes, in place of 16
#
# It checks the prestress (maxMag(D) falls with each pass and matches the
# reference) and the start of the FSI run (the solid inner-wall displacement,
# the fluid probe pressure, and FSI convergence within nOuterCorr). The
# prestressed wall starts in equilibrium with the initial 10 mmHg, so its
# displacement over the first steps must stay well below that of the
# unprestressed wall under the same pressure.
#
# The case needs PETSc and cartesianMesh. Like cerebralAneurysm, it is run by
# hand rather than from Alltest-regression.
# ============================================================

# ------------------------------------------------------------
# Regression settings
# ------------------------------------------------------------

N_PROCS=4
N_PASSES=3

# Regression end time for the copied case only: 5 time steps of 2.5e-4 s
REG_END_TIME=0.00125

# ------------------------------------------------------------
# Reference values and tolerances
# ------------------------------------------------------------

# Reference: the case's own answer on OpenFOAM-v2412 with 4 processes, which
# repeats to every digit written. The tolerances are wide enough to survive a
# compiler or PETSc version change and narrow enough to catch the material,
# the prestress or the coupling moving

# Prestress: maxMag(D) [m] of the runs that pass 1 and pass N_PASSES report,
# i.e. without prestress and with the Fpre of pass N_PASSES - 1
REF_PASS1_MAX_D=7.37231e-3
REF_PASSN_MAX_D=1.77080e-3
PRESTRESS_REL_TOL=1e-2

# FSI at REG_END_TIME: largest inner-wall displacement component [m] (2%) and
# the probe pressure [Pa], which falls by about 5 Pa from the initial 1333 Pa
# over these steps
REF_MAX_DISP=9.50797e-05
DISP_ABS_TOL=2e-6
REF_PROBE_P=1328.02
PROBE_P_ABS_TOL=1

# Blood density, to convert the kinematic probe pressure to Pa
RHO_FLUID=1060

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
HISTORY_FILE="prestress/prestressHistory"
FPRE_FILE="0/solid/Fpre"
DISP_FILE="postProcessing/0/solidDisplacementsinnerWall.dat"
FSI_FILE="postProcessing/0/fsiConvergenceData.dat"
PROBE_DIR="postProcessing/pressureProbe/fluid"

echo "============================================================"
echo "pulmonaryArteries prestressed FSI regression test"
echo "Processes                    = ${N_PROCS}"
echo "Prestress passes             = ${N_PASSES}"
echo "Regression end time          = ${REG_END_TIME}"
echo "Prestress maxMag(D) rel. tol = ${PRESTRESS_REL_TOL}"
echo "Max displacement abs. tol    = ${DISP_ABS_TOL}"
echo "Probe pressure abs. tol      = ${PROBE_P_ABS_TOL}"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    # Copy only the case inputs: results, logs and generated meshes from any
    # previous run of the tutorial are deliberately left behind
    for item in Allrun Allclean Allmesh 0 constant system prestress \
        plotInletVsPressure.gplt
    do
        cp -a "${SCRIPT_DIR}/${item}" "${CASE_DIR}/"
    done

    rm -f "${CASE_DIR}"/run.slurm "${CASE_DIR}"/prestress/run.slurm

    # Fewer prestress passes
    "${SOLIDS4FOAM_SED}" -i "s/^nPasses=[0-9]\+/nPasses=${N_PASSES}/" \
        "${CASE_DIR}/prestress/Allrun"

    # A few FSI time steps, written at the end so that reconstructPar has a
    # time to reconstruct
    "${SOLIDS4FOAM_SED}" -i \
        -e "s/^\(endTime[[:space:]]*\).*/\1${REG_END_TIME};/" \
        -e "s/^\(writeInterval[[:space:]]*\).*/\1${REG_END_TIME};/" \
        "${CASE_DIR}/system/controlDict"

    # Fewer processes, in every decomposeParDict
    for dict in \
        system/decomposeParDict \
        system/fluid/decomposeParDict \
        system/solid/decomposeParDict \
        prestress/system/decomposeParDict
    do
        "${SOLIDS4FOAM_SED}" -i \
            "s/^\(numberOfSubdomains[[:space:]]*\).*/\1${N_PROCS};/" \
            "${CASE_DIR}/${dict}"
    done
}

latest_numeric_time() {
    local file="$1"
    awk '
        ($1 + 0) == $1 { time = $1 }
        END {
            if (time != "") print time
        }
    ' "$file"
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
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun parallel > "${ALLRUN_LOGFILE}" 2>&1 ) \
        || true
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# The case requires PETSc and cartesianMesh: when either is unavailable Allrun
# writes a declared skip message, which is the only valid reason to skip the
# checks. Anything else that leaves the expected output missing or incomplete
# is a failure.
if solids4Foam::regressionCaseSkipped "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "Skipping regression checks because the tutorial skipped in this environment"
    exit 0
fi

# Allrun prints Done after the FSI run, the reconstruction and the plot
if ! grep -q "^Done" "${CASE_DIR}/${ALLRUN_LOGFILE}"; then
    echo "FAIL: Allrun did not complete:"
    tail -n 3 "${CASE_DIR}/${ALLRUN_LOGFILE}" | sed 's/^/      /'
    exit 1
fi

for file in "${HISTORY_FILE}" "${FPRE_FILE}" "${DISP_FILE}" "${FSI_FILE}"; do
    if [[ ! -f "${CASE_DIR}/${file}" ]]; then
        echo "FAIL: the case did not run in this environment:"
        echo "      ${file} is missing and the tutorial did not declare a skip"
        echo "      (see ${CASE_DIR}/${ALLRUN_LOGFILE})"
        exit 1
    fi
done

# The probe writes one directory per start time; the run starts from 0
PROBE_FILE=$(ls "${CASE_DIR}/${PROBE_DIR}"/*/p 2>/dev/null | head -1 || true)

if [[ -z "${PROBE_FILE}" ]]; then
    echo "FAIL: the probe pressure is missing (see ${CASE_DIR}/${ALLRUN_LOGFILE})"
    exit 1
fi

for file in "${CASE_DIR}/${DISP_FILE}" "${CASE_DIR}/${FSI_FILE}" "${PROBE_FILE}"
do
    t=$(latest_numeric_time "${file}" || true)

    if [[ -z "${t}" ]] \
    || ! awk "BEGIN {exit !(${t} + 0 >= ${REG_END_TIME} - 1e-9)}"
    then
        echo "FAIL: ${file#"${CASE_DIR}"/} stops at t = ${t:-none}, short of"
        echo "      the requested end time ${REG_END_TIME}: the case did not"
        echo "      complete"
        exit 1
    fi
done

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

# maxMag(D), the second column, of each pass in prestressHistory
prestress_max_d() {
    awk '!/^#/ && NF >= 2 {print $2}' "${CASE_DIR}/${HISTORY_FILE}"
}

extract_max_displacement() {
    # solidDisplacements writes:
    #     time minX minY minZ maxX maxY maxZ
    # The largest displacement component magnitude on the final line is used
    tail -1 "${CASE_DIR}/${DISP_FILE}" | \
    awk '
    {
        maxD = 0
        for (i = 2; i <= 7; i++)
        {
            d = ($i < 0 ? -$i : $i)
            if (d > maxD) maxD = d
        }
        print maxD
    }'
}

extract_max_fsi_correctors() {
    # fsiConvergenceData writes: time nFsiCorrectors
    awk '!/^#/ && NF >= 2 && $2 > n {n = $2} END {print n + 0}' \
        "${CASE_DIR}/${FSI_FILE}"
}

extract_probe_pressure() {
    # probes writes: time p, kinematic
    awk -v rho="${RHO_FLUID}" '!/^#/ && NF >= 2 {p = $2} END {print p*rho}' \
        "${PROBE_FILE}"
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

mapfile -t pass_max_d < <(prestress_max_d)
n_passes_run=${#pass_max_d[@]}

max_disp=$(extract_max_displacement)
max_fsi_corr=$(extract_max_fsi_correctors)
probe_p=$(extract_probe_pressure)
n_outer_corr=$(awk '$1 == "nOuterCorr" {gsub(";", "", $2); print $2}' \
    "${CASE_DIR}/constant/fsiProperties")

if [[ -z "${max_disp}" || -z "${max_fsi_corr}" || -z "${probe_p}" \
   || -z "${n_outer_corr}" ]]
then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

check_abs() {
    local name="$1" value="$2" ref="$3" tol="$4"
    local diff
    diff=$(abs "$(awk "BEGIN {print ${value} - (${ref})}")")

    if awk "BEGIN {exit !(${diff} <= ${tol})}"; then
        printf "PASS: %s = %.6g (ref %.6g, Δ = %.3g)\n" \
            "${name}" "${value}" "${ref}" "${diff}"
    else
        printf "FAIL: %s = %.6g (ref %.6g, Δ = %.3g > %.3g)\n" \
            "${name}" "${value}" "${ref}" "${diff}" "${tol}"
        failures=$((failures + 1))
    fi
}

check_rel() {
    local name="$1" value="$2" ref="$3" tol="$4"
    check_abs "${name}" "${value}" "${ref}" \
        "$(awk "BEGIN {r = ${ref}; print (r < 0 ? -r : r)*${tol}}")"
}

# Prestress: one line per pass, falling maxMag(D)
if (( n_passes_run == N_PASSES )); then
    echo "PASS: prestressHistory has ${n_passes_run} passes"
else
    echo "FAIL: prestressHistory has ${n_passes_run} passes, not ${N_PASSES}"
    failures=$((failures + 1))
fi

if (( n_passes_run >= 2 )); then
    falling=yes
    for (( i = 1; i < n_passes_run; i++ )); do
        if ! awk "BEGIN {exit !(${pass_max_d[i]} < ${pass_max_d[i-1]})}"; then
            falling=no
        fi
    done

    if [[ "${falling}" == yes ]]; then
        echo "PASS: prestress maxMag(D) falls with each pass:" \
             "$(printf '%.4g ' "${pass_max_d[@]}")"
    else
        echo "FAIL: prestress maxMag(D) does not fall with each pass:" \
             "$(printf '%.4g ' "${pass_max_d[@]}")"
        failures=$((failures + 1))
    fi
fi

if (( n_passes_run >= 1 )); then
    check_rel "prestress pass 1 maxMag(D)" "${pass_max_d[0]}" \
        "${REF_PASS1_MAX_D}" "${PRESTRESS_REL_TOL}"
    check_rel "prestress pass ${n_passes_run} maxMag(D)" \
        "${pass_max_d[n_passes_run-1]}" "${REF_PASSN_MAX_D}" \
        "${PRESTRESS_REL_TOL}"
fi

# FSI: every time step converged within nOuterCorr
if (( max_fsi_corr < n_outer_corr )); then
    echo "PASS: every FSI time step converged" \
         "(at most ${max_fsi_corr} of ${n_outer_corr} correctors)"
else
    echo "FAIL: an FSI time step used all ${n_outer_corr} correctors"
    failures=$((failures + 1))
fi

# FSI: the prestressed wall starts in equilibrium, so it moves far less than
# the unprestressed wall did under the same pressure in pass 1
if (( n_passes_run >= 1 )) \
&& awk "BEGIN {exit !(${max_disp} < 0.1*${pass_max_d[0]})}"
then
    printf "PASS: the prestressed wall stays near the imaged geometry"
    printf " (%.4g m, against %.4g m without prestress)\n" \
        "${max_disp}" "${pass_max_d[0]:-0}"
else
    printf "FAIL: the prestressed wall moves %.4g m, not well below the" \
        "${max_disp}"
    printf " %.4g m without prestress\n" "${pass_max_d[0]:-0}"
    failures=$((failures + 1))
fi

check_abs "max inner-wall displacement [m]" "${max_disp}" \
    "${REF_MAX_DISP}" "${DISP_ABS_TOL}"
check_abs "probe pressure [Pa]" "${probe_p}" \
    "${REF_PROBE_P}" "${PROBE_P_ABS_TOL}"

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
