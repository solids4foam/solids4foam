set terminal postscript enhanced color solid

set output "deflection.ps"
set xlabel "Time, t [s]"
set ylabel "Dx [m]"
set y2label "Dy [m]"
set grid

set y2tics

plot [0.2:] "./history/0/point.dat" using 1:2 axis x1y1 title "Dx" with lines, \
"./history/0/point.dat" using 1:3 axis x1y2 title "Dy" with lines

#set output
#set terminal x11
