# Plot the normal stress along one interface for each mid-layer thickness,
# against the digitised curves of Ioka et al. (2007) where they are supplied.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; variant='...'; thicknesses='0 0.002 0.005 0.01';
#                interface='13'; output='...'" plotThicknessStudy.gnuplot

# The terminal follows the extension of the output file: pdf, or png for the
# figures of the tutorial README
if (output[strlen(output) - 2:] eq "png") {
    set terminal pngcairo dashed enhanced font "Times,32" size 1920,1080
    labelFont = ",38"
    keyFont = ",28"
    pointSize = 3.0
    pointLineWidth = 2.5
    curveLineWidth = 5.0
} else {
    set terminal pdfcairo dashed enhanced size 5in,3.5in
    labelFont = ",12"
    keyFont = ",8"
    pointSize = 0.35
    pointLineWidth = 0.5
    curveLineWidth = 1.5
}
set output output

set grid
set logscale xy
set xrange [5e-5:1]
set yrange [0.6:3]
set format x "10^{%L}"
set ytics (0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 3)
set xlabel "r / W" font labelFont
set ylabel "{/Symbol s}_{xx} / {/Symbol s}_0, interface ".interface\
    font labelFont
set key inside right top opaque box font keyFont

hasData(file) = system("[ -s '".file."' ] && echo 1 || echo 0") + 0

plot_command = "plot "
separator = ""

do for [i = 1:words(thicknesses)] {
    thickness = word(thicknesses, i)

    # Without a mid-layer, the two interfaces are the single interface 12
    name = (thickness eq "0") ? "12" : interface
    file = sprintf("%s/thickness_%s_h%s_interface%s.txt", profiles, variant,\
        thickness, name)
    if (hasData(file)) {
        plot_command = plot_command . separator . sprintf("'%s' using 1:2 \
            every 3 with points pt %d ps %g lw %g lc %d title 'h/W = %s'",\
            file, 2*i + 2, pointSize, pointLineWidth, i, thickness)
        separator = ", "
    }

    curve = sprintf("%s/reference_h%s_interface%s.txt", profiles, thickness,\
        interface)
    if (hasData(curve)) {
        plot_command = plot_command . separator . sprintf("'%s' using 1:2 \
            with lines lw %g lc %d title 'Ioka et al., h/W = %s'",\
            curve, curveLineWidth, i, thickness)
        separator = ", "
    }
}

eval plot_command
