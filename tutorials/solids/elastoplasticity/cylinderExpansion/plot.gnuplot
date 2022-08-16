set terminal pdfcairo enhanced font "Times-Roman,9" linewidth 4 rounded

set style line 80 lt 1 lw 2
set style line 81 lt 0  lw 2 # dashed
set style line 81 lt rgb "#808080"  # grey
set grid back linestyle 81
set border 3 back linestyle 80
set xtics nomirror
set ytics nomirror
set style line 1 lt 1
set style line 2 lt 2

set output "plot.pdf"

set encoding iso_8859_2
set xlabel "Inner Boundary Radius (in mm)"
set ylabel "Inner Boundary Radial Stress {/Symbol s}_{/Times-Italic rr} (in MPa)"
set key top right;

set xrange [10:85]
#set yrange [0:0.5]

# Parameters for Analytical solution
a0 = 10;
r0 = 10;
b0 = 20;

plot \
    -(0.5/sqrt(3)) \
    *log( ((r0/a0)**2 + (x/a0)**2 - 1)/((b0/a0)**2 + (x/a0)**2 - 1) ) w l \
    t "Analytical", \
    "postProcessing/0/solidStressesinner.dat" \
    u (1000*$1*(0.075/15) + 10):(-1e-6*$5) every 1000::2 w lp \
    lt 1 lw 0.5 pt 6 ps 2 t "Numerical"