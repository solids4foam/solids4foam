#!/bin/bash
# Continue an interrupted 4x run from its latest complete field write.
source /tmp/claude-1000/-home-philipc-OpenFOAM-philipc-v2512-solids4foam-development--claude-worktrees-bridge-cse-011abXXq2tYCXUW1Vn7LZBYy/d2864dc3-2f71-41af-99cf-0cb70ce9f33b/scratchpad/ht_env.sh
cd /home/philipc/OpenFOAM/philipc-v2512/solids4foam.development/.claude/worktrees/bridge-cse_011abXXq2tYCXUW1Vn7LZBYy/tutorials/fluidSolidInteraction/HronTurek/verification/work/iqnils_mesh_4x || exit 1
echo "start $(date -Is)"
if grep -q -E '^End\s*$' log.solids4Foam; then echo "run already complete"; exit 0; fi
part=$(ls -d postProcessing.part* 2>/dev/null | wc -l); part=$((part+1))
cp -a postProcessing postProcessing.part$part
cp log.solids4Foam log.solids4Foam.part$part
# Drop a time directory that the interrupted job may have left half written
ref=$(ls processor0 | grep -E '^[0-9.]+$' | sort -g | tail -2 | head -1)
last=$(ls processor0 | grep -E '^[0-9.]+$' | sort -g | tail -1)
nref=$(find processor*/$ref -type f | wc -l); nlast=$(find processor*/$last -type f | wc -l)
if [ "$nlast" -lt "$nref" ]; then echo "removing incomplete $last ($nlast < $nref files)"; rm -rf processor*/$last; fi
# A total-strain law needs no incremental restart fields (the solver refuses
# to continue without them unless told so)
grep -q 'restart no;' constant/solid/solidProperties || python3 -c "import re;p='constant/solid/solidProperties';s=open(p).read();open(p,'w').write(re.sub(r'(Coeffs\"?\s*\{)', r'\g<1>\n    restart no;', s, count=1))"
echo "restarting from $(ls processor0 | grep -E '^[0-9.]+$' | sort -g | tail -1)"
mpirun -np 32 solids4Foam -parallel > log.solids4Foam 2>&1
echo "solver exit $? $(date -Is)"
