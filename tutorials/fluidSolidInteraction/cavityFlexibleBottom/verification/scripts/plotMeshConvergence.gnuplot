# Plot the computed steady displacement and interface force against the
# digitised Tukovic et al. (2018) mesh study.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; variant='...'; displacements='...';
#                forces='...'; output='...'" plotMeshConvergence.gnuplot

set terminal pdfcairo dashed enhanced size 5in,6in
set output output

data = sprintf("%s/%s_meshConvergence.txt", profiles, variant)

set datafile separator ","
set multiplot layout 2,1

set grid
set logscale x
set xlabel "Mesh spacing {/Symbol D}x (m)"
set key inside right bottom opaque box

set ylabel "Steady u_y at (4 -1 0.5) (m)"
plot displacements using 1:2 with linespoints pt 7 ps 0.6 lc rgb "#333333" \
         title "Tukovic et al. (2018)", \
     data using 1:2 with linespoints pt 6 ps 0.8 lc rgb "#d7191c" \
         title "solids4foam"

set ylabel "Steady interface F_y (N)"
plot forces using 1:($2*20) with linespoints pt 7 ps 0.6 lc rgb "#333333" \
         title "Tukovic et al. (2018)", \
     data using 1:3 with linespoints pt 6 ps 0.8 lc rgb "#d7191c" \
         title "solids4foam"

unset multiplot
