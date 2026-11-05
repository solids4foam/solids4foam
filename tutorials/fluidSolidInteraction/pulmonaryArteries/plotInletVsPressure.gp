# Inlet volumetric flow rate vs pressure-probe pressure
# Usage: gnuplot plotInletVsPressure.gp
#
# Flow:     constant/fluid/inletFlow/humanFlow (m^3/s, period 0.86 s, repeated)
# Pressure: postProcessing/pressureProbe/fluid/*/p (kinematic, m^2/s^2),
#           all restart time directories merged; where runs overlap the
#           later restart overrides the earlier run,
#           converted to mmHg with rho = 1060 kg/m^3

rho    = 1060.0
mmHg   = 133.322
period = 0.86
flowFile  = 'constant/fluid/inletFlow/humanFlow'
probeDir   = 'postProcessing/pressureProbe/fluid'
joined    = 'inletFlow_probeP.dat'

# Interpolate the (periodic) inlet flow table onto the probe times
system(sprintf("c=1e30; for d in $(ls %s | sort -gr); do \
  awk -v c=$c '!/^#/ && $1+0 <= c+0' %s/$d/p; c=$d; done | sort -g | awk -v T=%g ' \
  NR==FNR { if (match($0, /\\([ \\t]*[-0-9.eE+]+[ \\t]+[-0-9.eE+]+/)) { \
              split(substr($0, RSTART+1, RLENGTH-1), a); n++; t[n]=a[1]; q[n]=a[2] } next } \
  /^#/ { next } \
  { tm = $1 - T*int($1/T); \
    for (i=1; i<n; i++) if (tm>=t[i] && tm<=t[i+1]) break; \
    Q = q[i] + (q[i+1]-q[i])*(tm-t[i])/(t[i+1]-t[i]); \
    print $1, Q, $2 }' %s - > %s", probeDir, probeDir, period, flowFile, joined))

set terminal pngcairo size 1400,550 enhanced font 'Arial,12'
set output 'inletFlow_vs_probePressure.png'
set multiplot layout 1,2

# Left: both quantities against time
set title 'Inlet flow rate and probe pressure vs time'
set xlabel 'Time [s]'
set ylabel 'Inlet flow rate [mL/s]'
set y2label 'Probe pressure [mmHg]'
set ytics nomirror
set y2tics
set grid
set key top left
plot joined u 1:($2*1e6) w l lw 2 lc rgb '#1f77b4' t 'Inlet flow rate', \
     joined u 1:($3*rho/mmHg) axes x1y2 w l lw 2 lc rgb '#d62728' t 'Probe pressure'

# Right: pressure-flow loop (colour = time)
unset y2label; unset y2tics; set ytics mirror
set title 'Probe pressure vs inlet flow rate'
set xlabel 'Inlet flow rate [mL/s]'
set ylabel 'Probe pressure [mmHg]'
set cblabel 'Time [s]'
set palette rgb 33,13,10
plot joined u ($2*1e6):($3*rho/mmHg):1 w l lw 2 lc palette notitle

unset multiplot
