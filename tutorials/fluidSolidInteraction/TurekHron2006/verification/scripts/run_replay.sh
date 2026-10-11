#!/bin/bash
# Usage: run_replay.sh <case dir>; study binary (as the 1x-4x replays)
source ~/ht_prod/env_study.sh
cd "$1" || exit 1
n=$(grep nprocs energy_setup.txt | awk "{print \$2}")
log=log.solids4Foam; [ -f $log ] && mv $log log.solids4Foam.prev.$(date +%s)
echo "start $(date -Is) n=$n host $(hostname) bin $(which solids4Foam)" >> run.out
mpirun -np $n solids4Foam -parallel > $log 2>&1
echo "exit $? end $(date -Is) log $log" >> run.out
