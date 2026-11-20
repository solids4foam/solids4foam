# Plot the relative L2 errors of w and M against the mesh spacing for the
# mesh variants, with first- and second-order guide slopes.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; variants='hex distHex tri'; output='...'"
#           plotConvergence.gnuplot

set terminal pdfcairo dashed enhanced size 5in,3.5in
set output output

set grid
set logscale xy
set format x "%g"
set format y "%g"
set xlabel "Mesh spacing (m)"
set ylabel "Relative L2 error"
set key outside right top box

n = words(variants)
colours = "#d7191c #2c7bb6 #1a9641"

first = sprintf("%s/%s_errors.txt", profiles, word(variants, 1))
stats first using 1:2 nooutput
coarse_h = STATS_max_x
coarse_e = STATS_max_y
fine_h = STATS_min_x

set xrange [fine_h/1.6:coarse_h*1.6]
set yrange [coarse_e/400:coarse_e*40]
set xtics (coarse_h, coarse_h/2, coarse_h/4, coarse_h/8)

plot for [i=1:n] sprintf("%s/%s_errors.txt", profiles, word(variants, i)) \
         using 1:2 with linespoints pt 6 ps 0.8 lc rgb word(colours, i) \
         title sprintf("%s, w", word(variants, i)), \
     for [i=1:n] sprintf("%s/%s_errors.txt", profiles, word(variants, i)) \
         using 1:3 with linespoints dt 2 pt 8 ps 0.8 lc rgb word(colours, i) \
         title sprintf("%s, M", word(variants, i)), \
     coarse_e*(x/coarse_h)**1 with lines dt 2 lc rgb "#777777" \
         title "1st order", \
     coarse_e*(x/coarse_h)**2 with lines dt 3 lc rgb "#777777" \
         title "2nd order"
