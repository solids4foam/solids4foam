set terminal pdfcairo enhanced color solid

set output "displacement.pdf"

set xlabel "Time (in s)"
set ylabel "Vertical displacement of container apex (in m)"
set grid
set key top right
set ytics 0.4

plot "postProcessing/0/solidPointDisplacement_disp.dat" u 1:3 w l notitle
