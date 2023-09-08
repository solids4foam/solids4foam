set terminal pdfcairo enhanced color solid

set output "force.pdf"
set xlabel "Time, t [s]"
set ylabel "Fx [N]"
set y2label "Fy [N]"
set grid

set y2tics

plot [0.1:] \
    "< sed s/[\\(\\)]//g forces/0/forces.dat" using 1:($2+$5) axis x1y1 title "Fx" with lines, \
    "< sed s/[\\(\\)]//g forces/0/forces.dat" using 1:($3+$6) axis x1y2 title "Fy" with lines