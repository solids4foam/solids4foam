set terminal pdfcairo enhanced color solid

set output "displacement.pdf"
set xlabel "Time, t [s]"
set ylabel "Radial displacement [m]"
set grid

plot \
    "./postProcessing/0/solidPointDisplacement_pointDispX.dat" using 1:2 title "u_r at {/Symbol q} = 0" with lines lw 2, \
    "./postProcessing/0/solidPointDisplacement_pointDispY.dat" using 1:3 title "u_r at {/Symbol q} = 90 deg" with lines lw 2
