#!/bin/bash -l
#SBATCH --account=p201141
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --time=07:00:00
# Smooth graded mesh family (ht_mesh_smooth.py): fluid-only trajectory replays packed on one
# 128-core node, each case in its own job step with the ranks in its energy_setup.txt (the sum
# must not exceed 128). After each run: checkMesh of the deformed mesh at the last write.
# Optional converged-solve check: if <work>/<case>_conv exists, it is restarted from the
# second-last write of <case> (27.6 s) once <case> has finished.
# Usage: sbatch -J htsm -o <work>/slurm-%j.out run_smooth_replays_meluxina.slurm <work> <case> [<case> ...]
# Binary: the study binary of the ring+box 8x' replay (solids4foam aee0c8e3, OpenFOAM v2512
# source build, foss-2024a GCC 13.3).
source /project/home/p201287/philip/solids4foam-hronturek-v2512/env-meluxina.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
work=$1; shift

run_case() {
    cd "$work/$1" || return 1
    local n
    n=$(awk '$1 == "nprocs" {print $2}' energy_setup.txt)
    { echo start $(date -Is) host $(hostname) n=$n job $SLURM_JOB_ID; echo bin $(which solids4Foam)
      echo s4f $(git -C /project/home/p201287/philip/solids4foam-hronturek-v2512 rev-parse HEAD)
      gcc --version | head -1; } >> run.out
    srun --exact -n "$n" --cpu-bind=cores solids4Foam -parallel > log.solids4Foam 2>&1
    echo exit $? end $(date -Is) >> run.out
    local last
    last=$(\ls -1d processor0/[0-9]* | sed 's|processor0/||' | sort -g | tail -1)
    srun --exact -n "$n" --cpu-bind=cores checkMesh -parallel -time "$last" -noFunctionObjects \
        > log.checkMesh.$last 2>&1
}

for c in "$@"; do
    (
        run_case "$c"
        if [ -d "$work/${c}_conv" ]; then
            t=$(\ls -1d "$work/$c"/processor0/[0-9]* | sed 's|.*/||' | sort -g | tail -2 | head -1)
            for p in "$work/$c"/processor*; do
                \cp -r "$p/$t" "$work/${c}_conv/$(basename "$p")/"
            done
            run_case "${c}_conv"
        fi
    ) &
done
wait
