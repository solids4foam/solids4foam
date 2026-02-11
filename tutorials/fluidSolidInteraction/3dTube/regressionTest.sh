#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Beam-in-cross-flow FSI regression test
# ============================================================

# ------------------------------------------------------------
# Regression tolerances
# ------------------------------------------------------------

DISP_MAX_TOL=1e-5      # max displacement absolute tolerance
FORCE_MEAN_TOL=1e-3    # mean force tolerance

# Number of samples from end of force.dat to average
FORCE_AVG_SAMPLES=50

# Reference values
REF_MAX_DISP=0.000149739
REF_MEAN_FORCE=-0.0474451

# Log files
ALLRUN_LOGFILE="log.Allrun"

# Data files
DISP_FILE="postProcessing/0/solidPointDisplacement_displacement.dat"
FORCE_FILE="postProcessing/fluid/forces/0/force.dat"

echo "============================================================"
echo "Beam-in-cross-flow FSI regression test"
echo "Max displacement difference < ${DISP_MAX_TOL}"
echo "Mean force difference       < ${FORCE_MEAN_TOL}"
echo "============================================================"
echo

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
    ./Allclean > /dev/null 2>&1 || true
    ./Allrun > "${ALLRUN_LOGFILE}" 2>&1
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

# OpenFOAM variant compatibility
mkdir -p postProcessing/fluid/forces/0
(
    cd postProcessing/fluid/forces/0

    # foam-extend writes forces to a 'forces' sub-directory so we will create a
    # link
    if [[ ! -e force.dat && -f ../../../../forces/0/forces.dat ]]; then
        ln -s ../../../../forces/0/forces.dat force.dat
    fi

    # OpenFOAM.org creates forces.dat instead of force.dat
    if [[ ! -e force.dat && -f forces.dat ]]; then
        ln -s forces.dat force.dat
    fi
)

# ------------------------------------------------------------
# Extract helpers
# ------------------------------------------------------------

extract_max_displacement() {
    awk '{print $3}' "${DISP_FILE}" | sort -g | tail -1
}

extract_mean_force_tail() {
    tail -n "${FORCE_AVG_SAMPLES}" "${FORCE_FILE}" | \
    awk '
    {
        # Remove parentheses
        gsub(/[()]/, "", $0)

        # After cleanup, fields are:
        # $1 = time
        # $2 = Fx (both formats)
        sum += $2
        n++
    }
    END {
        if (n > 0) print sum/n
    }'
}

abs() {
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

# ------------------------------------------------------------
# Extract values
# ------------------------------------------------------------

max_disp=$(extract_max_displacement)
mean_force=$(extract_mean_force_tail)

if [[ -z "${max_disp}" || -z "${mean_force}" ]]; then
    echo "FAIL: Could not extract regression quantities"
    exit 1
fi

disp_diff=$(awk "BEGIN {print ${max_disp} - ${REF_MAX_DISP}}")
disp_diff_abs=$(abs "${disp_diff}")

force_diff=$(awk "BEGIN {print ${mean_force} - ${REF_MEAN_FORCE}}")
force_diff_abs=$(abs "${force_diff}")

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

if awk "BEGIN {exit !(${disp_diff_abs} < ${DISP_MAX_TOL})}"; then
    printf "PASS: max displacement = %.6g (Δ = %.3g)\n" \
        "${max_disp}" "${disp_diff_abs}"
else
    printf "FAIL: max displacement = %.6g (Δ = %.3g)\n" \
        "${max_disp}" "${disp_diff_abs}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${force_diff_abs} < ${FORCE_MEAN_TOL})}"; then
    printf "PASS: mean force = %.6g (Δ = %.3g)\n" \
        "${mean_force}" "${force_diff_abs}"
else
    printf "FAIL: mean force = %.6g (Δ = %.3g)\n" \
        "${mean_force}" "${force_diff_abs}"
    failures=$((failures + 1))
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
