set terminal postscript enhanced

set output "x-force.eps"
set xlabel "Time, s"
set ylabel "Lift force, N"
set grid

set key left bottom

set border lw 2

set size 0.8,0.65

plot [2:] "< sed s/[\\(\\)]//g ./forces/0/forces.dat" using 1:(2*($2+$5)/0.447575692) axis x1y1 notitle with lines lw 2 lt 15


