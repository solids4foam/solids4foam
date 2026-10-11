# Overlay the closing analysis window of a verification run on the Turek-Hron
# FSI3 reference history. Both histories are written by the verification
# driver with time measured from the last maximum of uy, so the phases align.
#
# Usage:
#   gnuplot -e "run='postProcessing/iqnils_mesh_2x_history.csv'; \
#               reference='postProcessing/reference_history.csv'; \
#               output='postProcessing/iqnils_mesh_2x_history.png'; \
#               label='iqnils_mesh_2x'" scripts/plotPeriodicHistory.gnuplot
#
# Optional: benchmark='FSI2' and xmin=-1.1 for the FSI2 benchmark.

if (!exists("benchmark")) benchmark = "FSI3"
if (!exists("xmin")) xmin = -0.6

set terminal pngcairo size 1400,1000 enhanced font "Helvetica,12"
set output output
set datafile separator ","
set grid
set key top left
set xrange [xmin:0.02]
set xlabel "Time from the last u_y maximum, t [s]"

set multiplot layout 2,2 title "Turek-Hron ".benchmark.": ".label." vs the reference history" noenhanced

set ylabel "u_x(A) [mm]"
plot reference using 1:($4*1000) with lines lw 2 lc rgb "black" title "Turek-Hron reference", \
     run using 1:($4*1000) with lines lw 2 lc rgb "red" title label noenhanced

set ylabel "u_y(A) [mm]"
plot reference using 1:($5*1000) with lines lw 2 lc rgb "black" title "Turek-Hron reference", \
     run using 1:($5*1000) with lines lw 2 lc rgb "red" title label noenhanced

set ylabel "Drag [N/m]"
plot reference using 1:2 with lines lw 2 lc rgb "black" title "Turek-Hron reference", \
     run using 1:2 with lines lw 2 lc rgb "red" title label noenhanced

set ylabel "Lift [N/m]"
plot reference using 1:3 with lines lw 2 lc rgb "black" title "Turek-Hron reference", \
     run using 1:3 with lines lw 2 lc rgb "red" title label noenhanced

unset multiplot
