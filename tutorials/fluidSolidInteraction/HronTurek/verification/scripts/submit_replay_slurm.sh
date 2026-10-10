#!/bin/bash
# Usage: submit.sh <case dir> [partition]; submits the replay case to Slurm on xenosim
case=$(readlink -f "$1"); part=${2:-main}
n=$(grep nprocs "$case/energy_setup.txt" | awk '{print $2}')
sbatch -p "$part" -n "$n" -N 1 -J "h2_$(basename "$case")" -o "$case/slurm.out" \
    --wrap "bash $HOME/ht_2nd/scripts/run_replay_variant.sh $case"
