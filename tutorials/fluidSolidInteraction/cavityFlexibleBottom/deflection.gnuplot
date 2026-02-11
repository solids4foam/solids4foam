set terminal pdfcairo enhanced color solid

set output "deflection.pdf"
set xlabel "Time, t [s]"
set ylabel "u_x [m]"
set y2label "u_y [m]"
set grid

set y2tics
set yrange[-0.0001:-0.0002]
set y2range[-0.3:-0.05]

plot [2:] \
    "./postProcessing/0/solidPointDisplacement_pointDisp.dat" using 1:2 axis x1y1 title "u_x" with lines lw 2, \
    "./postProcessing/0/solidPointDisplacement_pointDisp.dat" using 1:3 axis x1y2 title "u_y" with lines lw 2
