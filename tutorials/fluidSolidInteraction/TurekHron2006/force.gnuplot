# Only this case's own output (postProcessing, or forces on foam-extend):
# verification/work holds other runs
set terminal pdfcairo enhanced color solid

set output "force.pdf"
set xlabel "Time, t [s]"
set ylabel "Fx [N/m]"
set y2label "Fy [N/m]"
set grid

set y2tics

plot [0.01:] \
    "< sed s/[\\(\\)]//g `find postProcessing forces -name 'force.dat' 2>/dev/null`" using 1:($2)/0.015 axis x1y1 title "Fx" with lines, \
    "< sed s/[\\(\\)]//g `find postProcessing forces -name 'force.dat' 2>/dev/null`" using 1:($3)/0.015 axis x1y2 title "Fy" with lines