# Plot the relative error of the free-edge deflection against the mesh
# spacing, with first- and second-order guide slopes.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; variant='...'; output='...';
#               metric1='...'; metric2='...'" plotConvergence.gnuplot

set terminal pdfcairo dashed enhanced size 5in,3.5in
set output output

data = sprintf("%s/%s_errors.txt", profiles, variant)

set grid
set logscale xy
set format x "%g"
set format y "%g"
set ytics (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
set xlabel "Mesh spacing (m)"
set ylabel "Relative error of the deflection"
set key inside right bottom opaque box

# Anchor the guide slopes on the coarsest point, which is the largest
# spacing and the largest error.
stats data using 1:2 nooutput
coarse_h = STATS_max_x
coarse_e = STATS_max_y
fine_h = STATS_min_x

# Keep the axes on the data, so the guide slopes do not stretch the plot.
fine_e = STATS_min_y
set xrange [fine_h/1.6:coarse_h*1.6]
set yrange [fine_e/3:coarse_e*2]
set xtics (coarse_h, coarse_h/2, coarse_h/4, coarse_h/8)

plot data using 1:2 with linespoints pt 6 ps 0.8 lc rgb "#d7191c" \
         title metric1, \
     data using 1:3 with linespoints pt 8 ps 0.8 lc rgb "#2c7bb6" \
         title metric2, \
     coarse_e*(x/coarse_h)**1 with lines dt 2 lc rgb "#777777" \
         title "1st order", \
     coarse_e*(x/coarse_h)**2 with lines dt 3 lc rgb "#777777" \
         title "2nd order"
