#!/bin/bash
# Coupled FSI3 8x'/s16x (752 896 fluid cells, 1680 x 96 solid), dt 5e-4, strict IQN-ILS tol 1e-5,
# nOuterCorr 30, no reuse, cold start to 6 s, on one MeluXina node (128 ranks), as a chain of
# 48 h segments (writes every 0.1 s, purgeWrite 2, stopAt writeNow 40 min before each segment
# ends). The first segment builds the case. Extra segments exit at once if the run is complete.
#   bash submit_8x_meluxina.sh [segments, default 5]
cd /project/home/p201136/philip-ht || exit 1
C=c_8x_s16x_dt5e-4_tol1e-5; nseg=${1:-5}
[ -e $C ] && { echo "$C exists: chain continuation segments by hand"; exit 1; }
j=$(SETUP="8 16 128 --write-every 0.1 --purge-write 2" sbatch --parsable -N 1 -J ht8_s1 \
    -o logs/$C-%j.out run_coupled_meluxina.slurm $C) || exit 1
echo "segment 1: $j"
for s in $(seq 2 "$nseg"); do
    j=$(sbatch --parsable -N 1 -J ht8_s$s -d afterany:$j -o logs/$C-%j.out run_coupled_meluxina.slurm $C) || exit 1
    echo "segment $s: $j"
done
