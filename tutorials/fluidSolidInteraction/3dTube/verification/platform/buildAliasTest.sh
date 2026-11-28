#!/bin/bash
# Build and run the standalone tmp-aliasing reproducer with several flag sets.
# Source an OpenFOAM environment first; run from this directory.
set -u
INC="-I$FOAM_SRC/OpenFOAM/lnInclude -I$FOAM_SRC/OSspecific/POSIX/lnInclude"
BASE="-std=c++17 -m64 -DWM_DP -DWM_LABEL_SIZE=32 -DNoRepository -ftemplate-depth-100"
LIB="-L$FOAM_LIBBIN -lOpenFOAM -ldl -lm -Wl,-rpath,$FOAM_LIBBIN"
echo "OpenFOAM $WM_PROJECT_VERSION; $(g++ --version | head -1)"
for v in "O3:-O3" "O2:-O2" "O1:-O1" "O0:-O0" \
         "O3-norestrict:-O3 -D__restrict__=" "O3-nofma:-O3 -ffp-contract=off" \
         "O3-native:-O3 -march=native"
do
    name=${v%%:*}; flags=${v#*:}
    g++ $BASE $flags $INC Test-tmpDotAlias.C -o Test-tmpDotAlias_$name $LIB \
        || { echo "$name: build failed"; continue; }
    echo "$name ($flags): $(./Test-tmpDotAlias_$name | grep 'max|')"
    rm -f Test-tmpDotAlias_$name
done
