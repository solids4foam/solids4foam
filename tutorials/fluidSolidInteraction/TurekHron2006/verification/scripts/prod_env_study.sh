# Study binary of the replay campaign: solids4foam aee0c8e35 (src = 27c9bfdf7), v2512, GCC 11.4 -O3
source /etc/profile.d/petsc.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source ~/bin/load-openfoam v2512 >/dev/null
export FOAM_USER_LIBBIN=$HOME/ht_variants/ht_build/lib
export FOAM_USER_APPBIN=$HOME/ht_variants/ht_build/bin
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$(echo "$LD_LIBRARY_PATH" | tr ":" "\n" | grep -v "philipc-v2512/platforms" | paste -sd:)
export PATH=$FOAM_USER_APPBIN:$(echo "$PATH" | tr ":" "\n" | grep -v "philipc-v2512/platforms" | paste -sd:)
