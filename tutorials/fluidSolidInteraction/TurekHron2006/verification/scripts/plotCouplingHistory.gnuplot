# Overlay the IQN-ILS and Robin-Neumann transients after the coupling start.
#
# Usage:
#   gnuplot -e "data='postProcessing/coupling_comparison_1x_history.csv'; \
#               output='postProcessing/coupling_comparison_1x_history.png'" \
#           scripts/plotCouplingHistory.gnuplot

set terminal pngcairo size 1400,1000 enhanced font "Helvetica,12"
set output output
set datafile separator ","
set grid
set key top left
set xlabel "Time, t [s]"

set multiplot layout 2,2 title "Turek-Hron FSI3: IQN-ILS vs Robin-Neumann coupling"

set ylabel "u_x(A) [mm]"
plot data using 1:($2*1000) with lines lw 2 lc rgb "black" title "IQN-ILS", \
     data using 1:($6*1000) with lines lw 2 dt 2 lc rgb "red" title "Robin-Neumann"

set ylabel "u_y(A) [mm]"
plot data using 1:($3*1000) with lines lw 2 lc rgb "black" title "IQN-ILS", \
     data using 1:($7*1000) with lines lw 2 dt 2 lc rgb "red" title "Robin-Neumann"

set ylabel "Drag [N/m]"
plot data using 1:4 with lines lw 2 lc rgb "black" title "IQN-ILS", \
     data using 1:8 with lines lw 2 dt 2 lc rgb "red" title "Robin-Neumann"

set ylabel "Lift [N/m]"
plot data using 1:5 with lines lw 2 lc rgb "black" title "IQN-ILS", \
     data using 1:9 with lines lw 2 dt 2 lc rgb "red" title "Robin-Neumann"

unset multiplot
