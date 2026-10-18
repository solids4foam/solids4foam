#!/bin/bash
#------------------------------------------------------------------------------
# MPI scaling calibration: the same number of time steps, from the same
# developed state of a completed (decomposed) case, with each number of
# ranks in turn, within one SLURM allocation so that all the runs see the
# same node conditions.
#
# Usage (from the verification directory), e.g.
#     sbatch --nodes=1 --ntasks=8 --cpus-per-task=1 \
#         --mem=16G --time=02:00:00 --export=ALL \
#         --output=work/calibration/scaling-%j.out \
#         scripts/mpi_scaling.sh work/calibration/cutLink-L4-dt0p001 300 "1 2 4 8"
#
# Results: <base>/../scaling/ranks-<n>/ and <base>/../scaling/scaling.csv
#------------------------------------------------------------------------------
set -euo pipefail
export USER="${USER:-$(id -un)}"
export HOME="${HOME:-$(getent passwd "$USER" | cut -d: -f6)}"

BASE="$(cd "$1" && pwd)"
NSTEPS="$2"
RANKS="${3:-1 2 4 8}"
SKIP=50    # time steps excluded from the per-step wall time (start-up)
OUT="$(dirname "$BASE")/scaling"

# Optional environment set-up script (see README.md); otherwise the
# environment of the submitting shell is used (sbatch --export=ALL)
if [ -n "${S4F_VERIFY_ENV:-}" ]
then
    source "$S4F_VERIFY_ENV"
    if declare -F s4f_check_binaries > /dev/null
    then
        s4f_check_binaries
    fi
fi
command -v solids4Foam > /dev/null || { echo "solids4Foam not found" >&2; exit 1; }

echo "job=${SLURM_JOB_ID:-} host=$(hostname) ntasks=${SLURM_NTASKS:-?}" \
     "base=$BASE nsteps=$NSTEPS ranks=$RANKS"

# Developed state of the base case, reconstructed once
cd "$BASE"
if [ -d processor0 ]
then
    reconstructPar -latestTime > log.reconstructPar.scaling 2>&1
fi
T0="$(ls -1d [0-9]* | sort -g | tail -1)"
[ "$T0" != "0" ] || { echo "ERROR: no developed state in $BASE"; exit 1; }
DT="$(awk '$1 == "deltaT" {gsub(";", "", $2); print $2}' system/controlDict)"
T1="$(python3 -c "print(f'{$T0 + $NSTEPS*$DT:.10g}')")"
echo "start time $T0, deltaT $DT, end time $T1"

mkdir -p "$OUT"
CSV="$OUT/scaling.csv"
echo "ranks,nsteps,solver_wall_s,s_per_step_after_$SKIP,final_Fx,final_Fy" > "$CSV"

for n in $RANKS
do
    CASE="$OUT/ranks-$n"
    rm -rf "$CASE"
    mkdir -p "$CASE"
    cp -a "$BASE/constant" "$BASE/system" "$CASE/"
    cp -a "$BASE/$T0" "$CASE/"
    cd "$CASE"

    sed -i \
        -e "s/^startFrom .*;/startFrom       latestTime;/" \
        -e "s/^endTime .*;/endTime         $T1;/" \
        -e "s/^writeControl .*;/writeControl    timeStep;/" \
        -e "s/^writeInterval .*;/writeInterval   1000000;/" \
        system/controlDict

    cat > system/decomposeParDict <<EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}

numberOfSubdomains  $n;

method          simple;

simpleCoeffs
{
    n           ($n 1 1);
}
EOF

    if [ "$n" -gt 1 ]
    then
        decomposePar -force > log.decomposePar 2>&1
        RUN=(mpirun -np "$n" --bind-to none solids4Foam -parallel)
    else
        RUN=(solids4Foam)
    fi

    START=$(date +%s.%N)
    "${RUN[@]}" > log.solids4Foam 2>&1
    END=$(date +%s.%N)

    WALL=$(python3 -c "print(f'{$END - $START:.2f}')")
    STEPS=$(grep -c '^Time = ' log.solids4Foam)
    PER_STEP=$(grep '^ExecutionTime' log.solids4Foam | awk -v s=$SKIP \
        'NR == s {a = $7} END {if (NR > s) printf "%.4f", ($7 - a)/(NR - s)}')
    FORCE=$(grep 'Immersed body cylinder: force' log.solids4Foam | tail -1 \
        | sed -E 's/.*force \(([^ ]+) ([^ ]+) .*/\1,\2/')
    echo "$n,$STEPS,$WALL,$PER_STEP,$FORCE" | tee -a "$CSV"
    cd "$BASE"
done

echo "Results: $CSV"
