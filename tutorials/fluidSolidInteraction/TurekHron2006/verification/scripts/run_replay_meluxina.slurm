#!/bin/bash -l
#SBATCH --account=p201136
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
# Usage: sbatch -J <name> -o <case>/slurm.out run_replay.slurm <case>
# Study binary on MeluXina: solids4foam aee0c8e3 (same commit as the xenosim study binary),
# source-built OpenFOAM v2512, foss-2024a (GCC 13.3)
source /project/home/p201287/philip/solids4foam-hronturek-v2512/env-meluxina.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd $1 || exit 1
log=log.solids4Foam; [ -f $log ] && log=log.solids4Foam.$(date +%s)
{ echo start $(date -Is) host $(hostname) n=128; echo bin $(which solids4Foam); echo s4f $(git -C /project/home/p201287/philip/solids4foam-hronturek-v2512 rev-parse HEAD); gcc --version | head -1; lscpu | grep "Model name"; } >> run.out
srun -n 128 --cpu-bind=cores solids4Foam -parallel > $log 2>&1
echo exit $? end $(date -Is) log $log >> run.out
