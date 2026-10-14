#!/bin/bash -l
#SBATCH --account=p201136
#SBATCH --partition=cpu
#SBATCH --qos=short
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
# Timing of coupled FSI3 cases: runs each given case (ht_coupled.py, endTime set short) to its endTime
# one after the other (or all at once with CONCURRENT=1), with the rank count of its coupled_setup.txt.
#   sbatch -N <nodes> -J ht_timing -o <dir>/slurm-%j.out run_timing_meluxina.slurm <case> [<case> ...]
source /project/home/p201136/philip-ht/env-ht8.sh
cases=(); for c in "$@"; do cases+=("$(readlink -f "$c")"); done
for c in "${cases[@]}"; do
    cd "$c" || continue
    n=$(awk '$1=="cores"{print $2}' coupled_setup.txt)
    echo "start $(date -Is) job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES n=$n s4f $(git -C "$S4F_ROOT" rev-parse --short HEAD)" >> run.out
    if [ -n "$CONCURRENT" ]; then
        ( srun --exact -n "$n" --cpu-bind=cores solids4Foam -parallel > log.solids4Foam.$SLURM_JOB_ID 2>&1
          echo "exit $? end $(date -Is)" >> run.out ) &
    else
        srun -n "$n" --cpu-bind=cores solids4Foam -parallel > log.solids4Foam.$SLURM_JOB_ID 2>&1
        echo "exit $? end $(date -Is)" >> run.out
    fi
done
wait
