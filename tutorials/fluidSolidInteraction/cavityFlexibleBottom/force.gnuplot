set terminal pdfcairo enhanced color solid

set output "force.pdf"
set xlabel "Time, t [s]"
set ylabel "Fx [N]"
set y2label "Fy [N]"
set grid

set y2tics
set yrange[-100:0]
set y2range[-1e4:0]

plot [0.01:] \
    "< sed s/[\\(\\)]//g `find . -name 'force.dat'`" using 1:2 axis x1y1 title "Fx" with lines lw 2, \
    "< sed s/[\\(\\)]//g `find . -name 'force.dat'`" using 1:3 axis x1y2 title "Fy" with lines lw 2