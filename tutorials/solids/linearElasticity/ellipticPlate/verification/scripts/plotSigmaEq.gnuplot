# Plot the sampled equivalent stress profiles against the digitised
# Demirdzic et al. (1997) reference curve.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; reference='...'; variant='...';
#                levels='1 2 3 4'; output='...'" plotSigmaEq.gnuplot

set terminal pdfcairo dashed enhanced size 5in,3.5in
set output output

set grid
set xrange [0:90]
set yrange [1:4]
set xlabel "Angle from the x symmetry plane (degrees)"
set ylabel "Equivalent (von Mises) stress (MPa)"
set key inside right bottom opaque box

plot_command = sprintf("plot '%s' using 1:2 with lines lw 2 lc rgb '#333333' \
    title 'Demirdzic et al. (1997)'", reference)

do for [i = 1:words(levels)] {
    level = word(levels, i)
    file = sprintf("%s/%s_mesh%s_sigmaEq.txt", profiles, variant, level)
    plot_command = plot_command . sprintf(", '%s' using 1:($2/1e6) \
        with linespoints pt 6 ps %g title 'mesh %s'", \
        file, 1.1 - 0.15*i, level)
}

eval plot_command
