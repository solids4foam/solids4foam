set terminal pdfcairo enhanced color solid

set output "deflection.pdf"
set xlabel "Time, t [s]"
set ylabel "dx [m]"
set y2label "dy [m]"
set grid

set y2tics

plot [2:] \
    "./postProcessing/0/solidPointDisplacement_pointDisp.dat" using 1:2 axis x1y1 title "Ux" with lines, \
    "./postProcessing/0/solidPointDisplacement_pointDisp.dat" using 1:3 axis x1y2 title "Uy" with lines
