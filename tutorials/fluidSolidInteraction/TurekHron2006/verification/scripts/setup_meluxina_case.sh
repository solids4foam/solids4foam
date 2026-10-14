#!/bin/bash
# Build one coupled FSI3 case on MeluXina with ht_coupled.py (production defaults).
#   setup_meluxina_case.sh <out dir> <fluid level> <solid level> <ranks> [extra ht_coupled.py args]
source /project/home/p201136/philip-ht/env-ht8.sh
out=$1 lev=$2 slev=$3 n=$4; shift 4
V=$S4F_ROOT/tutorials/fluidSolidInteraction/TurekHron2006
python3 $V/verification/scripts/ht_coupled.py --tutorial $V --level "$lev" --solid-level "$slev" \
    --delta-t 5e-4 --end 6 --cores "$n" --out "$out" "$@"
