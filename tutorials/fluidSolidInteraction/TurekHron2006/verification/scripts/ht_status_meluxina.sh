#!/bin/bash
# htstat.sh <case>...: last time, steps, clock, wall per step and iterations per step over the last 50 steps
for c in "$@"; do
  l=$(ls -t $c/log.solids4Foam* 2>/dev/null | head -1); [ -z "$l" ] && { echo "$c: no log"; continue; }
  awk -v c="$c" '
    /^Time = [0-9.e+-]+$/ {t=$3}
    /^Time = .*iteration:/ {it++}
    /^ExecutionTime/ {k++; T[k]=$7; I[k]=it}
    /FOAM FATAL|Floating point exception| nan |-nan/ {err=$0}
    END{ m=(k>50)?k-50:1;
      printf "%s: t=%s steps=%d clock=%ds last%d: %.2f s/step, %.2f it/step %s\n", c, t, k, T[k], k-m, (k>m)?(T[k]-T[m])/(k-m):0, (k>m)?(I[k]-I[m])/(k-m):0, err }' "$l"
  tail -1 $c/run.out 2>/dev/null | cut -c1-90
done
