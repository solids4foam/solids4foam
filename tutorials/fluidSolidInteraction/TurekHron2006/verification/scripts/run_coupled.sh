#!/bin/bash
# Usage: run_coupled.sh <case dir>; production build (PR #554), v2512, GCC 11.4 -O3, xenosim
source ~/ht_prod/env_prod.sh
cd "$1" || exit 1
n=$(awk '$1=="cores"{print $2}' coupled_setup.txt)
{ echo "start $(date -Is) n=$n host $(hostname)"; echo "solids4Foam $(which solids4Foam)";
  echo "commit $(git -C ~/ht_prod/s4f rev-parse HEAD)"; echo "WM_OPTIONS $WM_OPTIONS $(g++ --version | head -1)"; } >> run.out
log=log.solids4Foam; [ -f $log ] && log=log.solids4Foam.$(date +%s)
if [ "$n" -gt 1 ]; then mpirun -np $n solids4Foam -parallel > $log 2>&1; else solids4Foam > $log 2>&1; fi
echo "exit $? end $(date -Is) log $log" >> run.out
