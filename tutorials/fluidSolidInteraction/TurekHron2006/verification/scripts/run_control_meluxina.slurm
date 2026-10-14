#!/bin/bash -l
#SBATCH --account=p201136
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
# Platform control and restart check on one MeluXina node (production build 769640289):
#   <dir>/ctl_1x  coupled 1x'/s2x, cold 0-6 s (rank count of the xenosim run)
#   <dir>/ctl_2x  coupled 2x'/s4x, cold 0-6 s (rank count of the xenosim run)
#   <dir>/rst_1x  coupled 1x'/s2x, cold 0-5 s, writes every 0.25 s kept (purgeWrite 0): the
#                 continuous run of the restart check
#   <dir>/rst_r4.5   restart of rst_1x from its 4.5 s write, run to 5.0 s (writes 4.75, 5.0)
#   <dir>/rst_r4.75  restart of rst_r4.5 from its 4.75 s write (a restart of a restart), to 5.0 s
# The restarts use the rank count of rst_1x and branch as soon as the parent has written.
#   sbatch -J ht_ctl -o <dir>/slurm-%j.out run_control_meluxina.slurm <dir> [extra case ...]
# Extra cases (ht_coupled.py dirs) run concurrently to their endTime.
source /project/home/p201136/philip-ht/env-ht8.sh
dir=$(readlink -f "$1"); shift
ranks() { awk '$1=="cores"{print $2}' "$1/coupled_setup.txt"; }
start() {  # start <case>: run in the background with the case's rank count
    local c=$1 n; n=$(ranks "$c")
    echo "start $(date -Is) job $SLURM_JOB_ID n=$n host $(hostname) s4f $(git -C "$S4F_ROOT" rev-parse HEAD) $(g++ --version | head -1)" >> "$c/run.out"
    ( cd "$c" && srun --exact -n "$n" --cpu-bind=cores solids4Foam -parallel > log.solids4Foam.$SLURM_JOB_ID 2>&1
      echo "exit $? end $(date -Is)" >> run.out ) &
}
waitfor() {  # waitfor <case> <time after the write>: the write is complete once the next step runs
    until grep -q "^Time = $2\$" "$1"/log.solids4Foam.$SLURM_JOB_ID 2>/dev/null; do
        grep -q "^exit" "$1/run.out" 2>/dev/null && return 1   # the parent stopped first
        sleep 60
    done
}
branch() {  # branch <parent> <time> <new case> <endTime>
    local p=$1 t=$2 c=$3 e=$4
    mkdir -p "$c"; cp -r "$p/constant" "$p/system" "$p/coupled_setup.txt" "$c/"
    for d in "$p"/processor*; do mkdir -p "$c/$(basename "$d")"; cp -r "$d/constant" "$d/$t" "$c/$(basename "$d")/"; done
    sed -i "s/^endTime .*/endTime         $e;/; s/^purgeWrite .*/purgeWrite      0;/" "$c/system/controlDict"
    echo "branched from $p at $t $(date -Is)" > "$c/run.out"
}
start "$dir/rst_1x"; start "$dir/ctl_1x"; start "$dir/ctl_2x"
for c in "$@"; do start "$(readlink -f "$c")"; done
if waitfor "$dir/rst_1x" 4.5005; then
    branch "$dir/rst_1x" 4.5 "$dir/rst_r4.5" 5; start "$dir/rst_r4.5"
    if waitfor "$dir/rst_r4.5" 4.7505; then
        branch "$dir/rst_r4.5" 4.75 "$dir/rst_r4.75" 5; start "$dir/rst_r4.75"
    fi
fi
wait
