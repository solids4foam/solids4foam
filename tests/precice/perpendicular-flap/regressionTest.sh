#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

#------------------------------------------------------------------------------
# License
#     This file is part of solids4foam, licensed under GNU General Public
#     License <http://www.gnu.org/licenses/>.
#
# Script
#     tests/precice/perpendicular-flap/regressionTest.sh
#
# Description
#     Regression checks for the preCICE perpendicular-flap tutorial with the
#     fluid-openfoam and solid-solids4foam participants.
#
#     Run from the tutorial directory by Alltest-precice, which has already run
#     the participants. The quantities checked come from the preCICE watch-point
#     on the flap tip, declared in precice-config.xml as
#
#         <watch-point mesh="Solid-Mesh" name="Flap-Tip" coordinate="0.0;1" />
#
#     and from the preCICE implicit-coupling iteration log. The flap-tip
#     x-displacement is the quantity used to compare the solid participants of
#     this tutorial against each other; see precice/tutorials#909.
#
# Reference values
#     Generated with tests/precice/generateReferences.sh in the CI container
#     (OpenFOAM-v2512), with the versions pinned in tests/precice/versions.txt
#     at the time of writing: preCICE 3.4.1, precice/tutorials 9afc7c84.
#     Regenerate after changing any pinned version.
#
# Tolerances
#     Repeat runs on one machine reproduce bit-identically, so the tolerances
#     are set by the spread across environments rather than by run-to-run
#     noise. The measured differences are:
#
#       OpenFOAM-v2512 (CI) vs OpenFOAM-v2606   ~2e-4 in both displacements,
#                                               4 in the iteration count
#       1% change in Young's modulus            2.9e-3 in the final and
#                                               1.5e-3 in the maximum
#
#     5e-4 therefore sits about twice above the version spread and about three
#     times below the signal from a small physical change, so the checks pass
#     when run locally against a different OpenFOAM version while still biting
#     on a real regression.
#
#------------------------------------------------------------------------------

CASE_DIR="$(pwd)"
SOLID_DIR="${CASE_DIR}/solid-solids4foam"

WATCHPOINT_FILE="${SOLID_DIR}/precice-Solid-watchpoint-Flap-Tip.log"

# ------------------------------------------------------------
# Reference values and tolerances
# ------------------------------------------------------------

# End time the truncated case is expected to reach, matching MAX_TIME in
# caseSetup
REG_END_TIME=1.0

# Flap-tip x-displacement at REG_END_TIME
REF_FINAL_DISP_X=0.0740134463
FINAL_DISP_X_TOL=5e-4

# Maximum |flap-tip x-displacement| over the run
REF_MAX_DISP_X=0.16607829
MAX_DISP_X_TOL=5e-4

# Total implicit coupling iterations over the run. This catches a convergence
# regression that leaves the displacement history close enough to pass the
# checks above.
REF_TOTAL_ITERATIONS=204
TOTAL_ITERATIONS_TOL=15

echo "============================================================"
echo "perpendicular-flap preCICE regression test"
echo "  Regression end time             = ${REG_END_TIME}"
echo "  Final x-displacement difference < ${FINAL_DISP_X_TOL}"
echo "  Max x-displacement difference   < ${MAX_DISP_X_TOL}"
echo "  Coupling iteration difference  <= ${TOTAL_ITERATIONS_TOL}"
echo "============================================================"
echo

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

abs()
{
    awk -v x="$1" 'BEGIN {print (x < 0 ? -x : x)}'
}

# Index of a named column in the watch-point header. The header is a comment
# line of the form
#     Time  Coordinate0  Coordinate1  Displacement0  ...
# Looking the name up rather than hard-coding an index keeps this working if
# the upstream watch-point gains or loses a data set.
columnIndex()
{
    local file="$1"
    local name="$2"

    awk -v name="${name}" '
        NR == 1 {
            for (i = 1; i <= NF; i++)
            {
                field = $i
                gsub(/^#/, "", field)
                if (field == name) { print i; exit }
            }
        }
    ' "${file}"
}

# ------------------------------------------------------------
# Check that the run produced usable output
# ------------------------------------------------------------

if [[ ! -f "${WATCHPOINT_FILE}" ]]
then
    echo "FAIL: watch-point file not found: ${WATCHPOINT_FILE}"
    exit 1
fi

DISP_X_COL=$(columnIndex "${WATCHPOINT_FILE}" "Displacement0")

if [[ -z "${DISP_X_COL}" ]]
then
    echo "FAIL: no Displacement0 column in ${WATCHPOINT_FILE}"
    echo "Header: $(head -1 "${WATCHPOINT_FILE}")"
    exit 1
fi

final_time=$(
    awk '($1 + 0) == $1 { time = $1 } END { if (time != "") print time }' \
        "${WATCHPOINT_FILE}"
)

if [[ -z "${final_time}" ]]
then
    echo "Skipping regression checks because the case produced no watch-point data"
    exit 0
fi

# A run that stopped early is reported as a skip rather than a failure, matching
# the convention of the tutorial regression tests. The driver reports the
# participant exit status separately.
if ! awk "BEGIN {exit !(${final_time} + 0 >= ${REG_END_TIME} - 1e-9)}"
then
    echo "Skipping regression checks because the case did not reach the"
    echo "requested end time (reached ${final_time}, expected ${REG_END_TIME})"
    exit 0
fi

# ------------------------------------------------------------
# Extract the quantities
# ------------------------------------------------------------

final_disp_x=$(
    awk -v col="${DISP_X_COL}" \
        '($1 + 0) == $1 { value = $col } END { print value }' \
        "${WATCHPOINT_FILE}"
)

max_disp_x=$(
    awk -v col="${DISP_X_COL}" '
        ($1 + 0) == $1 {
            value = $col < 0 ? -$col : $col
            if (value > maximum) maximum = value
        }
        END { print maximum + 0 }
    ' "${WATCHPOINT_FILE}"
)

# preCICE writes one line per time window, with the total iteration count in the
# "Iterations" column
ITERATIONS_FILE="${SOLID_DIR}/precice-Solid-iterations.log"
total_iterations=""

if [[ -f "${ITERATIONS_FILE}" ]]
then
    ITER_COL=$(columnIndex "${ITERATIONS_FILE}" "Iterations")

    if [[ -n "${ITER_COL}" ]]
    then
        total_iterations=$(
            awk -v col="${ITER_COL}" \
                '($1 + 0) == $1 { sum += $col } END { print sum + 0 }' \
                "${ITERATIONS_FILE}"
        )
    fi
fi

if [[ -z "${final_disp_x}" || -z "${max_disp_x}" ]]
then
    echo "FAIL: could not extract the regression quantities"
    exit 1
fi

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

failures=0

final_diff=$(abs "$(awk "BEGIN {print ${final_disp_x} - ${REF_FINAL_DISP_X}}")")
max_diff=$(abs "$(awk "BEGIN {print ${max_disp_x} - ${REF_MAX_DISP_X}}")")

if awk "BEGIN {exit !(${final_diff} < ${FINAL_DISP_X_TOL})}"
then
    printf "PASS: final flap-tip x-displacement = %.6g (ref %.6g, d = %.3g)\n" \
        "${final_disp_x}" "${REF_FINAL_DISP_X}" "${final_diff}"
else
    printf "FAIL: final flap-tip x-displacement = %.6g (ref %.6g, d = %.3g)\n" \
        "${final_disp_x}" "${REF_FINAL_DISP_X}" "${final_diff}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${max_diff} < ${MAX_DISP_X_TOL})}"
then
    printf "PASS: max flap-tip |x-displacement| = %.6g (ref %.6g, d = %.3g)\n" \
        "${max_disp_x}" "${REF_MAX_DISP_X}" "${max_diff}"
else
    printf "FAIL: max flap-tip |x-displacement| = %.6g (ref %.6g, d = %.3g)\n" \
        "${max_disp_x}" "${REF_MAX_DISP_X}" "${max_diff}"
    failures=$((failures + 1))
fi

if [[ -z "${total_iterations}" ]]
then
    echo "WARN: no coupling iteration data found, skipping that check"
else
    iter_diff=$(abs "$((total_iterations - REF_TOTAL_ITERATIONS))")

    if (( iter_diff <= TOTAL_ITERATIONS_TOL ))
    then
        printf "PASS: total coupling iterations = %d (ref %d, d = %d)\n" \
            "${total_iterations}" "${REF_TOTAL_ITERATIONS}" "${iter_diff}"
    else
        printf "FAIL: total coupling iterations = %d (ref %d, d = %d)\n" \
            "${total_iterations}" "${REF_TOTAL_ITERATIONS}" "${iter_diff}"
        failures=$((failures + 1))
    fi
fi

echo
if (( failures == 0 ))
then
    echo "============================================================"
    echo "perpendicular-flap preCICE regression test PASSED"
    echo "============================================================"
    exit 0
else
    echo "============================================================"
    echo "perpendicular-flap preCICE regression test FAILED (${failures} checks)"
    echo "============================================================"
    exit 1
fi
