# Plot the stress error norms against the mesh spacing, with first- and
# second-order guide slopes.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; variant='...'; output='...'"
#           plotErrorNorms.gnuplot

set terminal pdfcairo dashed enhanced size 5in,3.5in
set output output

data = sprintf("%s/%s_errors.txt", profiles, variant)

set grid
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set xlabel "Effective mesh spacing"
set ylabel "Relative stress error"
set key inside left top opaque box

# Anchor the guide slopes on the coarsest point, which is the largest
# spacing and the largest error.
stats data using 1:2 nooutput
coarse_h = STATS_max_x
coarse_e = STATS_max_y
fine_h = STATS_min_x

# Keep the axes on the data, so the guide slopes do not stretch the plot.
set xrange [fine_h/1.6:coarse_h*1.6]

plot data using 1:2 with linespoints pt 6 ps 0.8 lc rgb "#d7191c" \
         title "relative L_2", \
     data using 1:3 with linespoints pt 8 ps 0.8 lc rgb "#2c7bb6" \
         title "relative L_{/Symbol \245}", \
     coarse_e*(x/coarse_h)**1 with lines dt 2 lc rgb "#777777" \
         title "1st order", \
     coarse_e*(x/coarse_h)**2 with lines dt 3 lc rgb "#777777" \
         title "2nd order"
