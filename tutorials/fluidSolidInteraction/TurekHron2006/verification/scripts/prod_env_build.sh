# Production build: solids4foam verification/hronturek-production (incl. PR #554), v2512, GCC 11.4 -O3, xenosim
source /etc/profile.d/petsc.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source ~/bin/load-openfoam v2512 >/dev/null
unset FOAM_MODULE_LIBBIN FOAM_MODULE_APPBIN
export FOAM_USER_LIBBIN=$HOME/ht_prod/build/lib
export FOAM_USER_APPBIN=$HOME/ht_prod/build/bin
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$(echo "$LD_LIBRARY_PATH" | tr ":" "\n" | grep -v "philipc-v2512/platforms" | paste -sd:)
export PATH=$FOAM_USER_APPBIN:$(echo "$PATH" | tr ":" "\n" | grep -v "philipc-v2512/platforms" | paste -sd:)
