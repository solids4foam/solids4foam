set terminal pdfcairo enhanced color solid

set output "deflection.pdf"
set xlabel "Time, t [s]"
set ylabel "Centre point vertical deflection, u_y [m]"
set grid

#set yrange[-0.15:0.15]

plot "./postProcessing/0/solidPointDisplacement_pointDisp.dat" using 1:3 title "u_y" with lines lw 2
