# Plot one sampled set for each mesh level of the mesh study: the normal
# stress along an interface, or along the mid-line around the mid-layer.
#
# Invoked by the verification driver as:
#   gnuplot -e "profiles='...'; variant='...'; levels='1 2 3';
#                set='interface13'; thickness='0.01'; output='...'"
#           plotMeshStudy.gnuplot

# The terminal follows the extension of the output file: pdf, or png for the
# figures of the tutorial README
if (output[strlen(output) - 2:] eq "png") {
    set terminal pngcairo dashed enhanced font "Times,32" size 1920,1080
    labelFont = ",38"
    keyFont = ",28"
    pointSize = 2.0
    lineWidth = 3.5
    curveLineWidth = 5.0
} else {
    set terminal pdfcairo dashed enhanced size 5in,3.5in
    labelFont = ",12"
    keyFont = ",8"
    pointSize = 0.25
    lineWidth = 1.0
    curveLineWidth = 1.5
}
set output output

set grid
set key inside right top opaque box font keyFont

hasData(file) = system("[ -s '".file."' ] && echo 1 || echo 0") + 0

if (set eq "midLineZoom") {
    set xlabel "x / W" font labelFont
    set ylabel "{/Symbol s}_{xx} / {/Symbol s}_0, mid-line" font labelFont
    set format x "%.3f"
    set key inside left top
    style = "linespoints"
} else {
    set logscale xy
    set xrange [5e-5:1]
    set yrange [0.6:3]
    set format x "10^{%L}"
    set ytics (0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 3)
    set xlabel "r / W" font labelFont
    set ylabel "{/Symbol s}_{xx} / {/Symbol s}_0, interface ".set[10:11]\
        font labelFont
    style = "linespoints"
}

plot_command = "plot "
separator = ""

do for [i = 1:words(levels)] {
    level = word(levels, i)
    file = sprintf("%s/mesh_%s_mesh%s_%s.txt", profiles, variant, level, set)
    if (hasData(file)) {
        plot_command = plot_command . separator . sprintf("'%s' using 1:2 \
            with %s pt %d ps %g lw %g lc %d title 'mesh %s'",\
            file, style, 2*i + 2, pointSize, lineWidth, i, level)
        separator = ", "
    }
}

curve = sprintf("%s/reference_h%s_%s.txt", profiles, thickness, set)
if (hasData(curve)) {
    plot_command = plot_command . separator . sprintf("'%s' using 1:2 \
        with lines lw %g lc rgb '#333333' title 'Ioka et al. (2007)'",\
        curve, curveLineWidth)
}

eval plot_command
