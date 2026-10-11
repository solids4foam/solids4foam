#!/bin/bash
source ~/ht_prod/env_prod.sh
cd ~/ht_prod/s4f
echo "commit $(git rev-parse HEAD)"; echo "WM_OPTIONS $WM_OPTIONS"; g++ --version | head -1
./Allwmake -j 16 > ~/ht_prod/logs/build_prod.log 2>&1; echo "Allwmake exit $?"
ls -la $FOAM_USER_LIBBIN $FOAM_USER_APPBIN/solids4Foam
