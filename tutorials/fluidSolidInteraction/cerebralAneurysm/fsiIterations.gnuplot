# Plot the number of fluid-solid interaction (FSI) iterations per time-step.
#
# postProcessing/fsiResiduals.dat contains one line per FSI iteration, holding
# the time, the FSI outer-corrector index and the FSI residual. The awk filter
# below reduces this to the final corrector index reached at each time, i.e. the
# number of FSI iterations performed in that time-step.

fsiIters = "< awk 'NR > 1 { if (haveTime && $1 != time) print time, iters; time = $1; iters = $2; haveTime = 1 } END { if (haveTime) print time, iters }' postProcessing/fsiResiduals.dat"

set title "cerebralAneurysm: number of FSI iterations per time-step"
set xlabel "Time, t [s]"
set ylabel "Number of FSI iterations"
set grid
set key top right
set yrange [0:*]

set terminal pdfcairo enhanced color solid
set output "fsiIterations.pdf"
plot fsiIters u 1:2 w l lw 1 lt 1 lc "red" t "FSI iterations"

set terminal pngcairo enhanced color solid size 1000,600 font ",12"
set output "fsiIterations.png"
replot
