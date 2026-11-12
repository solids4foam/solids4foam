set terminal pdfcairo enhanced color solid

set output "deflection.pdf"
set xlabel "Time, t [s]"
set ylabel "Tail-tip displacement [m]"
set grid

plot \
    "./postProcessing/0/solidPointDisplacement_pointDisp.dat" using 1:3 title "Dy" with lines
