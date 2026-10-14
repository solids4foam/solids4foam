#!/bin/bash
# Coupled FSI3 8x'/s16x (752 896 fluid cells, 1680 x 96 solid), dt 5e-4, strict IQN-ILS tol 1e-5,
# nOuterCorr 30, no reuse, cold start to 6 s, solid PC hypre BoomerAMG (tolerances unchanged;
# equivalence checked on 1x', see README), one MeluXina node (128 ranks), account p201141,
# as an afterany chain of segments (writes every 0.1 s, purgeWrite 2, stopAt writeNow 40 min
# before each segment ends). The first segment builds the case; a segment that finds the run
# complete exits at once.
#   bash submit_8x_meluxina.sh long [segments, default 2]      6-day segments (QOS long, 1 running job)
#   bash submit_8x_meluxina.sh default [segments, default 4]   48 h segments
cd /project/home/p201136/philip-ht || exit 1
C=c_8x_s16x_dt5e-4_tol1e-5_hypre; qos=${1:?long or default}
case $qos in long) wall=6-00:00:00; nseg=${2:-2};; default) wall=48:00:00; nseg=${2:-4};; *) exit 1;; esac
[ -e $C ] && { echo "$C exists: chain continuation segments by hand"; exit 1; }
opts="-A p201141 -q $qos -t $wall -N 1"
j=$(SETUP="8 16 128 --write-every 0.1 --purge-write 2 --solid-pc hypre" sbatch --parsable $opts -J ht8_s1 \
    -o logs/$C-%j.out run_coupled_meluxina.slurm $C) || exit 1
echo "segment 1: $j"
for s in $(seq 2 "$nseg"); do
    j=$(sbatch --parsable $opts -J ht8_s$s -d afterany:$j -o logs/$C-%j.out run_coupled_meluxina.slurm $C) || exit 1
    echo "segment $s: $j"
done
