# Plot the computed reaction force histories against the digitised
# Neto et al. (2016) reference curves.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; referenceX='...'; referenceY='...';
#                variant='mu0.3'; levels='1 2 3'; displacement=31.5;
#                output='...'" plotReactionForces.gnuplot

set terminal pdfcairo dashed enhanced size 5in,6in
set output output

set multiplot layout 2,1

set grid
set xrange [0:1]
set xlabel "Normalised displacement"
set key inside left bottom opaque box

do for [component = 1:2] {
    if (component == 1) {
        set ylabel "Reaction force in x (N)"
        reference = referenceX
    } else {
        set ylabel "Reaction force in y (N)"
        reference = referenceY
    }

    plot_command = sprintf("plot '%s' using ($1/%g):2 with points pt 6 \
        lc rgb '#333333' title 'Neto et al. (2016)'", reference, displacement)

    do for [i = 1:words(levels)] {
        level = word(levels, i)
        file = sprintf("%s/%s_mesh%s_reaction.txt", profiles, variant, level)
        plot_command = plot_command . sprintf(", '%s' using ($1/%g):%d \
            with lines lw %g title 'mesh %s'", \
            file, displacement, component + 1, 2.4 - 0.4*i, level)
    }

    eval plot_command
}

unset multiplot
