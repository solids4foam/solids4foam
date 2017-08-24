set terminal postscript enhanced color solid

set output "force.ps"
set xlabel "Time, t [s]"
set ylabel "Fx [m]"
set y2label "Fy [m]"
set grid

set y2tics

plot [2:] "./history/0/force.dat" using 1:(-$2) axis x1y1 title "Ux" with lines, \
"./history/0/force.dat" using 1:(-$3) axis x1y2 title "Uy" with lines

#set output
#set terminal x11
