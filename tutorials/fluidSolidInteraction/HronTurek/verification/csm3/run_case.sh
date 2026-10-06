#!/bin/bash
# Run one CSM3 case: ./run_case.sh <name> <nx> <ny> <deltaT> <endTime> [overlayDir]
# Copies ./case to ./work/<name>, sets the mesh/time controls, runs blockMesh
# and solids4Foam. An optional overlay directory is copied over the case
# (e.g. to change fvSchemes or solidProperties for a sensitivity test).
set -e
name=$1; nx=$2; ny=$3; dt=$4; endTime=$5; overlay=$6
here=$(cd "$(dirname "$0")" && pwd)
dir=$here/work/$name
rm -rf "$dir"; mkdir -p "$here/work"; cp -r "$here/case" "$dir"
if [[ -n "$overlay" ]]; then cp -r "$overlay"/. "$dir"/; fi
cd "$dir"
sed -i "s/(105 6 1)/($nx $ny 1)/" system/blockMeshDict
sed -i -e "s/^deltaT .*/deltaT $dt;/" -e "s/^endTime .*/endTime $endTime;/" system/controlDict
blockMesh > log.blockMesh 2>&1
start=$(date +%s)
solids4Foam > log.solids4Foam 2>&1
echo "$(( $(date +%s) - start ))" > runtime.s
