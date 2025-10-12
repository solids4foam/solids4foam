set terminal pdfcairo enhanced color solid

set output "fsiConvergence.pdf"

set title "3dTube: Number of FSI iterations"
set xlabel "Time, t [s]"
set ylabel "Number of FSI iterations"
set grid
set key top right

set yrange [0:80]

plot \
    "postProcessing/0/fsiConvergenceData.dat" u 1:2 w lp ps 0.3 pt 5 lc "red" lw 0.5 lt 1 t "Number of iterations"
