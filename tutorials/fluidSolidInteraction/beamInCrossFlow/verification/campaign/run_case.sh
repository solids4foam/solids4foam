#!/bin/bash
# Mesh, decompose and run one case prepared by scripts/richter_case.py.
#
#   run_case.sh <case> <ranks> [launcher...]
#
# The launcher defaults to "mpirun -np <ranks>"; batch scripts pass e.g.
# "srun --exact -n 32 --cpu-bind=cores".  A case that already has processor
# directories is continued from its latest time (startFrom latestTime).
set -o pipefail
case_dir=$(cd "$1" && pwd) || exit 1
ranks=$2
shift 2
if [ $# -gt 0 ]; then
    launcher=("$@")
else
    launcher=(mpirun -np "$ranks")
fi
cd "$case_dir" || exit 1
if [ ! -d processor0 ] && [ ! -d processors"$ranks" ]; then
    blockMesh -region solid > log.blockMesh.solid 2>&1 || exit 1
    blockMesh -region fluid > log.blockMesh.fluid 2>&1 || exit 1
    checkMesh -region fluid > log.checkMesh.fluid 2>&1
    checkMesh -region solid > log.checkMesh.solid 2>&1
    if [ "$ranks" -gt 1 ]; then
        decomposePar -region fluid > log.decomposePar.fluid 2>&1 || exit 1
        decomposePar -region solid > log.decomposePar.solid 2>&1 || exit 1
    fi
    log=log.solids4Foam
else
    sed -i.bak -e 's/^startFrom .*/startFrom       latestTime;/' system/controlDict
    n=1
    while [ -e log.solids4Foam.restart$n ]; do n=$((n+1)); done
    log=log.solids4Foam.restart$n
fi
echo "$(date -Is) start $(hostname) ranks=$ranks log=$log" >> campaign_runs.txt
if [ "$ranks" -gt 1 ]; then
    "${launcher[@]}" solids4Foam -parallel > "$log" 2>&1
else
    solids4Foam > "$log" 2>&1
fi
status=$?
echo "$(date -Is) end status=$status" >> campaign_runs.txt
exit $status
