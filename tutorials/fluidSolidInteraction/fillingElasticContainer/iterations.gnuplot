set terminal pdfcairo enhanced color solid

set output "iterations.pdf"

set xlabel "Time (in s)"
set ylabel "Number of FSI iterations"
set grid
set key top right

plot "postProcessing/fsiResiduals.dat" u 1:2 w l lw 0.5 lt 1 notitle
