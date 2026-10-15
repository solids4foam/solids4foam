#!/bin/bash
# Health gate of a coupled FSI3 run (ht_coupled.py case): prints "OK ..." or "ABORT <reason>".
#   ht_health_gate.sh <case dir> [log file]
# Rules (calibrated on the healthy 4x' run at dt 5e-4: max Co <= 3.84 on all 10 000 steps, the
# first coupled steps converge in 13 iterations; and on the unstable 8x' run at dt 5e-4: max Co
# 8-34 from 0.03 s, the first coupled steps capped at residuals 1e-2):
#   1. Courant: after t = CO_T0 (0.05 s), abort if at least CO_N (20) of the last 100 steps have
#      max Co > CO_MAX (5).
#   2. Coupling start: abort if more than CAP_N (3) of the first 20 coupled steps reach nOuterCorr.
#   3. Coupling floor: abort if the median final residual of the last 20 coupled steps > RES_MAX (1e-4).
c=${1:?case}; log=${2:-$(ls -t "$c"/log.solids4Foam* 2>/dev/null | head -1)}
CO_MAX=${CO_MAX:-5} CO_N=${CO_N:-20} CO_T0=${CO_T0:-0.05} CAP_N=${CAP_N:-3} RES_MAX=${RES_MAX:-1e-4}
nout=$(awk '$1=="nOuterCorr"{gsub(";","",$2); print $2}' "$c/constant/fsiProperties")
tc=$(awk '$1=="couplingStartTime"{gsub(";","",$2); print $2}' "$c/constant/fsiProperties")
[ -n "$log" ] && [ -f "$log" ] || { echo "OK no log yet"; exit 0; }
co=$(tail -n 300000 "$log" | awk -v t0="$CO_T0" -v cm="$CO_MAX" '
    /^Time = [0-9.eE+-]+$/ {t=$3}
    /^Courant Number mean/ { if (t!=lt && t+0 > t0) { n++; v[n]=$6; tt[n]=t; lt=t } }
    END { k=0; mx=0; for (i=(n>100?n-99:1); i<=n; i++) { if (v[i]>cm) k++; if (v[i]>mx) mx=v[i] }
          printf "%d %g %s %d\n", k, mx, tt[n], n }')
read -r nco comax tlast nsteps <<< "$co"
if [ "${nsteps:-0}" -gt 0 ] && [ "$nco" -ge "$CO_N" ]; then
    echo "ABORT Courant: $nco of the last 100 steps up to t=$tlast have max Co > $CO_MAX (largest $comax)"; exit 1
fi
pp=$c/postProcessing
fsi=$( { for f in $(ls "$pp"/fsiResiduals.dat.upto_* 2>/dev/null | sort -t_ -k2 -g) "$pp/fsiResiduals.dat"; do
          [ -f "$f" ] && tail -n +2 "$f"; done; } | awk -v tc="$tc" -v no="$nout" -v cn="$CAP_N" -v rm="$RES_MAX" '
    $1+0 > tc+1e-9 { if (!($1 in it)) ord[++n]=$1; it[$1]=$2; r[$1]=$3 }
    END { n--;                                   # the newest step may be in progress
          cap=0; for (i=1; i<=n && i<=20; i++) if (it[ord[i]]+0 >= no) cap++
          m=0; for (i=(n>20?n-19:1); i<=n; i++) a[++m]=r[ord[i]]+0
          asort(a); med=(m? a[int((m+1)/2)] : 0)
          printf "%d %d %g\n", n, cap, med }')
read -r ncpl ncap rmed <<< "$fsi"
if [ "${ncpl:-0}" -ge 1 ] && [ "$ncap" -gt "$CAP_N" ]; then
    echo "ABORT coupling start: $ncap of the first $([ $ncpl -lt 20 ] && echo $ncpl || echo 20) coupled steps reached nOuterCorr=$nout"; exit 1
fi
if [ "${ncpl:-0}" -ge 20 ] && awk -v x="$rmed" -v y="$RES_MAX" 'BEGIN{exit !(x>y)}'; then
    echo "ABORT coupling floor: median final residual of the last 20 coupled steps $rmed > $RES_MAX"; exit 1
fi
echo "OK t=$tlast steps>Co$CO_MAX(last100)=$nco maxCo=$comax coupled=$ncpl capped(first20)=$ncap medianRes(last20)=$rmed"
