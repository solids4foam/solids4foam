# Plot the error norms against the analytical plate-with-hole solution as a
# function of the mesh spacing, with first- and second-order guide slopes.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; variant='...'; output='...'" \
#       plotConvergence.gnuplot
#
# The data file written by the driver has the columns
#   spacing_m displacement_l2_m displacement_linf_m stress_xx_l2_pa
#   stress_xx_linf_pa

set terminal pdfcairo dashed enhanced size 5in,3.5in
set output output

data = sprintf("%s/%s_convergence.txt", profiles, variant)

set grid
# The mesh family halves the spacing at every level, so a base-two log
# axis puts one tick on each level.
set logscale x 2
set logscale y
set xlabel "Effective cell spacing (m)"
set ylabel "Error norm, normalised by the coarsest mesh"
set key outside right top box
set format x "%.3f"
set format y "10^{%T}"

# Normalise every norm by its value on the coarsest mesh so that the
# displacement and stress errors share one set of guide slopes.
stats data using 1 nooutput
coarse_h = STATS_max
set xrange [0.7*STATS_min : 1.4*STATS_max]
stats data using (($1 == coarse_h) ? $2 : 1/0) nooutput
d_l2 = STATS_max
stats data using (($1 == coarse_h) ? $3 : 1/0) nooutput
d_linf = STATS_max
stats data using (($1 == coarse_h) ? $4 : 1/0) nooutput
s_l2 = STATS_max
stats data using (($1 == coarse_h) ? $5 : 1/0) nooutput
s_linf = STATS_max

first(h) = h/coarse_h
second(h) = (h/coarse_h)**2

plot data using 1:($2/d_l2) with linespoints pt 7 ps 0.5 lw 2 \
         title "D, L2", \
     data using 1:($3/d_linf) with linespoints pt 5 ps 0.5 lw 2 \
         title "D, LInf", \
     data using 1:($4/s_l2) with linespoints pt 9 ps 0.5 lw 2 \
         title "{/Symbol s}_{xx}, L2", \
     data using 1:($5/s_linf) with linespoints pt 11 ps 0.5 lw 2 \
         title "{/Symbol s}_{xx}, LInf", \
     data using 1:(first($1)) with lines dt 2 lc rgb "#777777" \
         title "first order", \
     data using 1:(second($1)) with lines dt 3 lc rgb "#333333" \
         title "second order"
