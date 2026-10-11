#!/bin/bash
# Usage: sbatch --ntasks=N run_cfd3.sh --level L --delta-t DT --end-time T --cores N [--tag X]
source /tmp/claude-1000/-home-philipc-OpenFOAM-philipc-v2512-solids4foam-development--claude-worktrees-bridge-cse-011abXXq2tYCXUW1Vn7LZBYy/d2864dc3-2f71-41af-99cf-0cb70ce9f33b/scratchpad/ht_env.sh
cd /home/philipc/OpenFOAM/philipc-v2512/solids4foam.development/.claude/worktrees/bridge-cse_011abXXq2tYCXUW1Vn7LZBYy/tutorials/fluidSolidInteraction/HronTurek/verification
echo "host $(hostname) start $(date -Is) args $*"
python3 scripts/hron_turek_cfd3.py "$@"
echo "exit $? end $(date -Is)"
