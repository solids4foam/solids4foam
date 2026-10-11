set terminal pdfcairo enhanced color solid

set output "deflection.pdf"
set xlabel "Time, t [s]"
set ylabel "Midpoint vertical displacement, u_y [m]"
set grid
set key bottom right

plot \
    "./postProcessing/0/solidPointDisplacement_midpoint.dat" using 1:3 \
    title "solids4foam" with lines lw 2
