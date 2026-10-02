set terminal pdfcairo enhanced color solid

set output "deflection.pdf"
set xlabel "Time, t [s]"
set ylabel "Roof centre vertical displacement [m]"
set grid
set key bottom right

plot [0:12] \
    "./postProcessing/0/solidPointDisplacement_pointDisp.dat" using 1:3 title "solids4foam" with lines lw 2, \
    "./reference/vonScheven2009_dz.dat" using 1:2 title "von Scheven (2009)" with lines lw 2 dt 2
