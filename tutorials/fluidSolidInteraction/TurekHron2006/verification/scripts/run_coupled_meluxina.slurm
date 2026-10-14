#!/bin/bash -l
#SBATCH --account=p201136
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
# One segment of a coupled FSI3 run (ht_coupled.py case) on MeluXina, production build
# 769640289 (src = 30d6e22f3 + PR #562), OpenFOAM v2512 source build, GCC 13.3.
#   sbatch -N <nodes> -J <name> -o <case>/slurm-%j.out run_coupled_meluxina.slurm <case> [margin min]
# Chain segments with --dependency=afterany:<previous job>. Each segment continues from the
# latest write (startFrom latestTime), and <margin> minutes (default 40) before the job's
# end time sets stopAt writeNow, so that the segment ends on a complete write
# (restart yes: D_0_0, D_0_0_0, DD_0_0, A and U are written, PR #562). A segment that
# finds the case at endTime exits at once.
# With SETUP="<fluid level> <solid level> <ranks> [ht_coupled.py args]" a missing case is first
# built with setup_meluxina_case.sh (modules are only available inside jobs).
source /project/home/p201136/philip-ht/env-ht8.sh
case=$(readlink -f "$1"); margin=${2:-40}
if [ ! -d "$case" ] && [ -n "$SETUP" ]; then
    mkdir -p "$(dirname "$case")"; cd "$(dirname "$case")"
    bash /project/home/p201136/philip-ht/setup_meluxina_case.sh "$case" $SETUP > "setup_$(basename "$case").log" 2>&1 || exit 1
fi
cd "$case" || exit 1
n=$(awk '$1=="cores"{print $2}' coupled_setup.txt)
end=$(awk '$1=="endTime"{gsub(";","",$2); print $2}' system/controlDict)
last=$(ls processor0 | grep -E '^[0-9][0-9.e+-]*$' | sort -g | tail -1)
if awk -v l="$last" -v e="$end" 'BEGIN{exit !(l >= e - 1e-9)}'; then
    echo "job $SLURM_JOB_ID: already at $last (endTime $end), nothing to do $(date -Is)" >> run.out; exit 0
fi
sed -i 's/^stopAt .*/stopAt          endTime;/' system/controlDict
# fsiResiduals.dat is rewritten from the start on every run: keep the previous segments'
[ -f postProcessing/fsiResiduals.dat ] && mv postProcessing/fsiResiduals.dat "postProcessing/fsiResiduals.dat.upto_$last"
log=log.solids4Foam.$SLURM_JOB_ID
{ echo "start $(date -Is) job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ($SLURM_JOB_NODELIST) n=$n from $last"
  echo "bin $(which solids4Foam)"; echo "s4f $(git -C "$S4F_ROOT" rev-parse HEAD)"
  echo "$WM_OPTIONS $(g++ --version | head -1)"; lscpu | grep "Model name"; } >> run.out
srun -n "$n" --cpu-bind=cores solids4Foam -parallel > "$log" 2>&1 &
pid=$!
jobend=$(date -d "$(scontrol show job "$SLURM_JOB_ID" | grep -oP 'EndTime=\K\S+')" +%s)
stopat=$((jobend - margin * 60)); stopped=
while kill -0 $pid 2>/dev/null; do
    sleep 60
    if [ -z "$stopped" ] && [ "$(date +%s)" -ge "$stopat" ]; then
        sed -i 's/^stopAt .*/stopAt          writeNow;/' system/controlDict
        stopped=1; echo "stopAt writeNow set $(date -Is)" >> run.out
    fi
done
wait $pid; rc=$?
sed -i 's/^stopAt .*/stopAt          endTime;/' system/controlDict
echo "exit $rc end $(date -Is) log $log last $(ls processor0 | grep -E '^[0-9][0-9.e+-]*$' | sort -g | tail -1)" >> run.out
