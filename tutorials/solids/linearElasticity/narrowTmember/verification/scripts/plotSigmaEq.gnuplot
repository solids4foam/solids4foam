# Plot the sampled equivalent stress profiles against the digitised
# Demirdzic et al. (1997) reference curve.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; reference='...'; variant='...';
#                levels='1 2 3 4'; output='...'" plotSigmaEq.gnuplot

set terminal pdfcairo dashed enhanced size 5in,3.5in
set output output

set grid
set xrange [0:270]
set xtics 45
set yrange [0:8]
set xlabel "Angle from (0, -1.5R) (degrees)"
set ylabel "Equivalent (von Mises) stress (MPa)"
set key inside left top opaque box

plot_command = sprintf("plot '%s' using 1:2 with linespoints pt 7 ps 0.5 \
    lc rgb '#333333' title 'Demirdzic et al. (1997)'", reference)

do for [i = 1:words(levels)] {
    level = word(levels, i)
    file = sprintf("%s/%s_mesh%s_sigmaEq.txt", profiles, variant, level)
    plot_command = plot_command . sprintf(", '%s' using 1:($2/1e6) \
        with lines lw %g title 'mesh %s'", file, 2.0 - 0.25*i, level)
}

eval plot_command
