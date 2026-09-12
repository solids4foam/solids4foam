# Plot the deformed mid-wall lines sampled from each mesh level.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; variant='...'; levels='1 2 3';
#                output='...'" plotMidLine.gnuplot

set terminal pdfcairo dashed enhanced size 4in,5in
set output output

set grid
set size ratio -1
set xrange [-14:0]
set yrange [-28:5]
set xlabel "x (mm)"
set ylabel "z (mm)"
set key outside center bottom

plot_command = "plot "

do for [i = 1:words(levels)] {
    level = word(levels, i)
    file = sprintf("%s/%s_mesh%s_midLine.txt", profiles, variant, level)
    separator = (i == 1) ? "" : ", "
    plot_command = plot_command . sprintf("%s'%s' using ($1*1e3):($3*1e3) \
        with linespoints pt 6 ps %g title 'mesh %s'", \
        separator, file, 1.1 - 0.15*i, level)
}

eval plot_command
