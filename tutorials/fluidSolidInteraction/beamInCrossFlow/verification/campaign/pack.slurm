#!/bin/bash -l
# Run several prepared cases concurrently on one or more MeluXina nodes.
#
#   sbatch -J <name> -N <nodes> -t <hh:mm:ss> pack.slurm <ranks> <case> [<case> ...]
#
# Every case gets <ranks> MPI ranks (srun --exact), so the sum over the cases
# must not exceed 128 x nodes.  Each case is meshed and decomposed on first
# use and continued from its latest time otherwise (see run_case.sh).
#SBATCH --account=p201287
#SBATCH --qos=default
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=128
#SBATCH --output=/project/home/p201287/philip/beam-richter-opus/jobs/%x-%j.out

BR=/project/home/p201287/philip/beam-richter-opus
source $BR/env.sh
ranks=$1
shift
echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NODELIST start $(date -Is)"
echo "solids4Foam: $(which solids4Foam)"
ldd $(which solids4Foam) | grep -E "solids4Foam|petsc" || true
echo "case commit: $(git -C $BR/solids4foam rev-parse HEAD)"
pids=()
for case_dir in "$@"; do
    $BR/solids4foam/tutorials/fluidSolidInteraction/beamInCrossFlow/verification/campaign/run_case.sh \
        "$case_dir" "$ranks" srun --exact -n "$ranks" --cpu-bind=cores &
    pids+=($!)
done
status=0
for i in "${!pids[@]}"; do
    wait "${pids[$i]}" || { echo "case ${@:$((i+1)):1} failed"; status=1; }
done
echo "end $(date -Is) status $status"
exit $status
