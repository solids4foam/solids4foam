# HT 8x' production build on MeluXina: solids4foam verification/hronturek-production 769640289
# (src = 30d6e22f3 + PR #562 / 541f5b008), source-built OpenFOAM v2512, foss-2024a (GCC 13.3)
source /project/home/p201287/philip/openfoam-v2512/env-modules.sh
source /project/home/p201287/philip/openfoam-v2512/OpenFOAM-v2512/etc/bashrc
S4F_ROOT=/project/home/p201136/philip-ht/s4f-8x
export PETSC_DIR="$EBROOTPETSC" PETSC_ARCH=
export EIGEN_DIR="$EBROOTEIGEN/include"
export FOAM_USER_LIBBIN="$S4F_ROOT/platforms/lib" FOAM_USER_APPBIN="$S4F_ROOT/platforms/bin"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH" PATH="$FOAM_USER_APPBIN:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
