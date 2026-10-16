# Midpoint displacement histories of a verification study against the
# published curves. Called by the verification driver with:
#   histories: CSV written by the driver (time, one column per member)
#   ncurves:   number of members
#   refdir:    verification/reference
#   output:    PNG file name

set terminal pngcairo enhanced size 1400,900 font ",12"
set output output
set datafile separator ","
set datafile commentschars "#"
set key outside right top
set grid
set xlabel "Time, t (s)"
set ylabel "Midpoint vertical displacement, u_y (m)"
set xrange [0:70]
set yrange [-0.02:0.32]

plot \
    refdir."/Valdes2007_fig5.12.csv" using 1:2 \
        with lines lw 3 lc rgb "#d62728" title "Valdes (2007)", \
    refdir."/KratosZorrilla.csv" using 1:2 \
        with lines lw 3 lc rgb "#17becf" title "Kratos (Zorrilla)", \
    refdir."/Mok2001_bild6.6.csv" using 1:2 \
        with lines lw 1 dt 2 lc rgb "#7f7f7f" title "Mok (2001), context", \
    for [i=2:ncurves+1] histories using 1:i \
        with lines lw 1.5 title columnheader(i)
