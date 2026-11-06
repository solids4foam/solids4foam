# Plot the errors against the exact solution for the three configurations
# of the annular plate, with first- and second-order guide slopes.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; variants='a b c'; output='...'"
#           plotConvergence.gnuplot

set terminal pdfcairo dashed enhanced size 5in,3.5in
set output output

set grid
set logscale xy
set format x "%g"
set format y "%g"
set ytics (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
set xlabel "Radial mesh spacing (m)"
set ylabel "Relative error"
set key inside right bottom opaque box

first = sprintf("%s/%s_errors.txt", profiles, word(variants, 1))
stats first using 1:2 nooutput
coarse_h = STATS_max_x
coarse_e = STATS_max_y
fine_h = STATS_min_x
fine_e = STATS_min_y
set xrange [fine_h/1.6:coarse_h*1.6]
set yrange [fine_e/3:coarse_e*2]
set xtics (coarse_h, coarse_h/2, coarse_h/4, coarse_h/8)

plot for [i = 1:words(variants)] \
         sprintf("%s/%s_errors.txt", profiles, word(variants, i)) \
         using 1:2 with linespoints pt 4+i ps 0.8 lw 1.5 \
         title sprintf("%s: free-edge w", word(variants, i)), \
     for [i = 1:words(variants)] \
         sprintf("%s/%s_errors.txt", profiles, word(variants, i)) \
         using 1:3 with linespoints pt 4+i ps 0.8 dt 4 \
         title sprintf("%s: L_2 of w", word(variants, i)), \
     coarse_e*(x/coarse_h)**1 with lines dt 2 lc rgb "#777777" \
         title "1st order", \
     coarse_e*(x/coarse_h)**2 with lines dt 3 lc rgb "#777777" \
         title "2nd order"
