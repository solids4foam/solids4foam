set terminal pdfcairo enhanced color solid

set output "deflection.pdf"
set xlabel "Time, t [s]"
set ylabel "Dx [m]"
set y2label "Dy [m]"
set grid

set y2tics

plot [0.2:] \
    "./postProcessing/0/solidPointDisplacement_displacement.dat" using 1:2 axis x1y1 title "Dx" with lines, \
    "./postProcessing/0/solidPointDisplacement_displacement.dat" using 1:3 axis x1y2 title "Dy" with lines