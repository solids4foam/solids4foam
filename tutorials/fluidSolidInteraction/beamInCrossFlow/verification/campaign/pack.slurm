#!/bin/bash -l
# Run several prepared cases concurrently on one or more SLURM nodes.
#
#   sbatch -A <account> -J <name> -N <nodes> -t <hh:mm:ss> \
#       -o <log dir>/%x-%j.out pack.slurm <ranks> <case> [<case> ...]
#
# Every case gets <ranks> MPI ranks (srun --exact), so the sum over the cases
# must not exceed the cores per node x nodes (128 on the MeluXina CPU nodes
# used for the campaign).  Each case is meshed and decomposed on first use and
# continued from its latest time otherwise (see run_case.sh).
#
# Environment:
#   S4F_ENV   optional script sourced first to set up OpenFOAM and solids4foam
#             (e.g. module loads and the OpenFOAM bashrc)
# The account, partition and log location are site specific and are passed
# to sbatch; the defaults below are only the generic resource layout.
#SBATCH --qos=default
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=128

if [ -n "${S4F_ENV:-}" ]; then
    # shellcheck disable=SC1090
    source "$S4F_ENV"
fi
# sbatch runs a copy of this script, so locate run_case.sh through the
# submission directory if this script's own directory does not hold it
script_dir=$(cd "$(dirname "$0")" && pwd)
[ -x "$script_dir/run_case.sh" ] || script_dir=${SLURM_SUBMIT_DIR:-$script_dir}
run_case="$script_dir/run_case.sh"
[ -x "$run_case" ] || { echo "run_case.sh not found; submit from campaign/"; exit 1; }
ranks=$1
shift
cases=("$@")
echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NODELIST start $(date -Is)"
echo "solids4Foam: $(which solids4Foam)"
ldd "$(which solids4Foam)" | grep -E "solids4Foam|petsc" || true
echo "case commit: $(git -C "$script_dir" rev-parse HEAD 2>/dev/null || echo unknown)"
pids=()
for case_dir in "${cases[@]}"; do
    "$run_case" \
        "$case_dir" "$ranks" srun --exact -n "$ranks" --cpu-bind=cores &
    pids+=($!)
done
status=0
for i in "${!pids[@]}"; do
    wait "${pids[$i]}" || { echo "case ${cases[$i]} failed"; status=1; }
done
echo "end $(date -Is) status $status"
exit $status
