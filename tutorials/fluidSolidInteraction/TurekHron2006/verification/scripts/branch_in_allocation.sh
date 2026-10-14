#!/bin/bash
# Branch a coupled case from one of its writes and run the branch as a job step inside a running
# single-node allocation (idle cores of the control job), optionally with the solid PETSc
# preconditioner replaced (equivalence check of a solver-path change at unchanged tolerances).
#   branch_in_allocation.sh <jobid> <parent case> <time> <new case> <endTime> [hypre]
jobid=$1 p=$(readlink -f "$2") t=$3 c=$4 e=$5 pc=$6
mkdir -p "$c"; c=$(readlink -f "$c")
cp -r "$p/constant" "$p/system" "$p/coupled_setup.txt" "$c/"
for d in "$p"/processor*; do mkdir -p "$c/$(basename "$d")"; cp -r "$d/constant" "$d/$t" "$c/$(basename "$d")/"; done
sed -i "s/^endTime .*/endTime         $e;/; s/^purgeWrite .*/purgeWrite      0;/" "$c/system/controlDict"
if [ "$pc" = hypre ]; then
    /usr/bin/python3.12 - "$c/system/solid/fvSolution" <<'PY'
import sys
p = sys.argv[1]; t = open(p).read()
old = "            pc_type bjacobi;\n            sub_pc_type lu;\n"
new = "".join(f"            {l}\n" for l in [
    "pc_type hypre;", "pc_hypre_type boomeramg;", 'pc_hypre_boomeramg_max_iter "1";',
    'pc_hypre_boomeramg_strong_threshold "0.7";', 'pc_hypre_boomeramg_grid_sweeps_up "1";',
    'pc_hypre_boomeramg_grid_sweeps_down "1";', 'pc_hypre_boomeramg_agg_nl "1";',
    'pc_hypre_boomeramg_agg_num_paths "1";', 'pc_hypre_boomeramg_max_levels "25";',
    "pc_hypre_boomeramg_coarsen_type HMIS;", "pc_hypre_boomeramg_interp_type ext+i;",
    'pc_hypre_boomeramg_P_max "1";', 'pc_hypre_boomeramg_truncfactor "0.3";'])
assert t.count(old) == 1
open(p, "w").write(t.replace(old, new))
PY
fi
n=$(awk '$1=="cores"{print $2}' "$c/coupled_setup.txt")
echo "branched from $p at $t (pc ${pc:-unchanged}) $(date -Is); step of job $jobid n=$n" > "$c/run.out"
cd "$c" && nohup srun --jobid="$jobid" --exact -n "$n" --cpu-bind=cores bash -lc \
    "source /project/home/p201136/philip-ht/env-ht8.sh; exec solids4Foam -parallel" \
    > log.solids4Foam.step 2>&1 &
