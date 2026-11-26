#!/bin/bash
#
# Derive the manufactured solution, generate the C++ header and check it.
#
# Usage: ./run_checks.sh [outDir] [rhoF] [rhoS] [full|flat] [muF]
#
#     ./run_checks.sh                     full variant, as used by the case
#     ./run_checks.sh generated_flat 1 1 flat
#
# The derivation runs 11 symbolic and 40-digit checks, writes
# <outDir>/mmsFsiExact.H and <outDir>/reference_values.txt, and
# check_generated.C then checks the generated header independently. The
# full variant with the default parameters reproduces ../src/mmsFsiExact.H
# and reference_values.txt; copy them there to update the case.
#
# Requires Python 3 with sympy and symengine (PYTHON, default python3) and a
# C++17 compiler (CXX, default c++). Each variant takes about six minutes.

cd "${0%/*}" || exit 1

out=${1:-generated}
rhoF=${2:-1}
rhoS=${3:-1}
variant=${4:-full}
muF=${5:-1/100}
PY=${PYTHON:-python3}

if ! "$PY" -c "import sympy, symengine" > /dev/null 2>&1
then
    echo "sympy and symengine are required: pip install sympy symengine"
    exit 1
fi

if ! "$PY" -u mms_fsi.py --outdir "$out" --rho-f "$rhoF" --rho-s "$rhoS" \
    --variant "$variant" --mu-f "$muF" > "log.$out.derive" 2>&1
then
    tail -5 "log.$out.derive"
    echo "DERIVATION CHECKS FAILED"
    exit 1
fi
grep -E "PASS|FAIL|checks passed" "log.$out.derive"

${CXX:-c++} -O1 -std=c++17 -I "$out" \
    -DREFFILE=\""$out"/reference_values.txt\" \
    check_generated.C -o "check_$out" \
 && "./check_$out" | tee "log.$out.check"
