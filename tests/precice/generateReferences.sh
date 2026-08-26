#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

#------------------------------------------------------------------------------
# License
#     This file is part of solids4foam, licensed under GNU General Public
#     License <http://www.gnu.org/licenses/>.
#
# Script
#     tests/precice/generateReferences.sh
#
# Description
#     Run the registered preCICE cases and print the quantities their
#     regressionTest.sh scripts compare against, so that the reference values
#     can be updated after a deliberate change to versions.txt or to the
#     upstream case configuration.
#
#     This does not edit the regressionTest.sh scripts: read the printed values,
#     run this twice to see the run-to-run spread, and set the tolerances by
#     hand with headroom above that spread.
#
# Usage
#     ./tests/precice/generateReferences.sh [CASE...]
#
#------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORK_DIR="${ROOT_DIR}/tutorialsTest-precice"
TUTORIALS_DIR="${WORK_DIR}/tutorials"

# Run the cases. Reference values are meaningless if the checks did not run, but
# a failing comparison is expected here: that is the whole point.
"${SCRIPT_DIR}/Alltest-precice" "$@" || true

echo
echo "============================================================"
echo "Reference values"
echo "============================================================"

for setup in "${SCRIPT_DIR}"/*/caseSetup
do
    [[ -f "${setup}" ]] || continue

    caseName="$(basename "$(dirname "${setup}")")"

    if (( $# > 0 )) && [[ ! " $* " == *" ${caseName} "* ]]
    then
        continue
    fi

    # shellcheck source=/dev/null
    ( source "${setup}"

      solidDir="${TUTORIALS_DIR}/${TUTORIAL}/solid-solids4foam"
      watchpoint=$(find "${solidDir}" -maxdepth 1 -name 'precice-*-watchpoint-*.log' \
                       2>/dev/null | head -1)
      iterations=$(find "${solidDir}" -maxdepth 1 -name 'precice-*-iterations.log' \
                       2>/dev/null | head -1)

      echo
      echo "${caseName} (max-time ${MAX_TIME:-full})"

      if [[ -z "${watchpoint}" ]]
      then
          echo "  no watch-point output found"
          exit 0
      fi

      awk '
          NR == 1 {
              for (i = 1; i <= NF; i++)
              {
                  field = $i; gsub(/^#/, "", field)
                  if (field == "Displacement0") col = i
              }
              next
          }
          ($1 + 0) == $1 {
              final = $col
              value = $col < 0 ? -$col : $col
              if (value > maximum) maximum = value
              time = $1
          }
          END {
              printf "  end time                  = %s\n", time
              printf "  REF_FINAL_DISP_X          = %.9g\n", final
              printf "  REF_MAX_DISP_X            = %.9g\n", maximum + 0
          }
      ' "${watchpoint}"

      if [[ -n "${iterations}" ]]
      then
          awk '
              NR == 1 {
                  for (i = 1; i <= NF; i++)
                  {
                      field = $i; gsub(/^#/, "", field)
                      if (field == "Iterations") col = i
                  }
                  next
              }
              ($1 + 0) == $1 { sum += $col }
              END { printf "  REF_TOTAL_ITERATIONS      = %d\n", sum + 0 }
          ' "${iterations}"
      fi
    )
done

echo
echo "Generated with:"
echo "  solids4foam  $(cd "${ROOT_DIR}" && git describe --always --dirty 2>/dev/null || echo unknown)"
echo "  OpenFOAM     ${WM_PROJECT_VERSION:-unknown}"
echo "  preCICE      $(pkg-config --modversion libprecice 2>/dev/null || echo unknown)"
echo "  tutorials    $(cd "${TUTORIALS_DIR}" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo
