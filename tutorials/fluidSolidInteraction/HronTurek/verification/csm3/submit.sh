#!/bin/bash
# sbatch wrapper: ./submit.sh <name> <nx> <ny> <deltaT> <endTime> [overlay]
# CSM3_ENV is a script that loads the OpenFOAM/solids4foam environment.
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$here/work"
sbatch --partition=main --ntasks=1 --job-name="csm3_$1" -o "$here/work/slurm_$1.log" \
  --wrap "bash -c 'source ${CSM3_ENV:-$HOME/bin/load-openfoam v2512} >/dev/null; which solids4Foam; $here/run_case.sh $*'"
