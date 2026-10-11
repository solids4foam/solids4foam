#!/bin/bash
# 4x mesh, interface tolerance 1e-6, restarted from the t = 4 s state of the
# 1e-5 run (fields with old-time levels), to bound the coupling error
source /tmp/claude-1000/-home-philipc-OpenFOAM-philipc-v2512-solids4foam-development--claude-worktrees-bridge-cse-011abXXq2tYCXUW1Vn7LZBYy/d2864dc3-2f71-41af-99cf-0cb70ce9f33b/scratchpad/ht_env.sh
cd /home/philipc/OpenFOAM/philipc-v2512/solids4foam.development/.claude/worktrees/bridge-cse_011abXXq2tYCXUW1Vn7LZBYy/tutorials/fluidSolidInteraction/HronTurek/verification/work/iqnils_mesh_4x_tol1e-06 || exit 1
echo "start $(date -Is) solids4Foam $(which solids4Foam)"
mpirun -np 32 solids4Foam -parallel > log.solids4Foam 2>&1
echo "solver exit $? $(date -Is)"
