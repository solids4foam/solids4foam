#!/bin/bash
# Usage: run_case.sh <case dir>; study binary (GCC 11.4 -O3, v2512) + libhtVariants
source /etc/profile.d/petsc.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source ~/bin/load-openfoam v2512 >/dev/null
export FOAM_USER_LIBBIN=$HOME/ht_variants/ht_build/lib
export FOAM_USER_APPBIN=$HOME/ht_variants/ht_build/bin
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$(echo "$LD_LIBRARY_PATH" | tr ":" "\n" | grep -v "philipc-v2512/platforms" | paste -sd:)
export PATH=$FOAM_USER_APPBIN:$(echo "$PATH" | tr ":" "\n" | grep -v "philipc-v2512/platforms" | paste -sd:)
cd "$1" || exit 1
n=$(grep nprocs energy_setup.txt | awk "{print \$2}")
echo "start $(date -Is) n=$n" > run.out
if [ "$n" -gt 1 ]; then mpirun -np $n solids4Foam -parallel > log.solids4Foam 2>&1; else solids4Foam > log.solids4Foam 2>&1; fi
echo "exit $? end $(date -Is)" >> run.out
