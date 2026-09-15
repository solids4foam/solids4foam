#!/usr/bin/gnuplot
#
# Compares the solids4foam vertical force and twisting torque on the
# sphere-displacement patch with the digitised reference data.
# The indentation phase runs from t = 0 to 4 and the twisting phase from
# t = 4 to 13 (20 degrees per unit time), so angle = 20*(t - 4).

set term pngcairo dashed font "Arial,18" size 1024,768

set linestyle 1 lt 6 lc rgb "black" lw 2 ps 1.5
set linestyle 2 lt 6 lc rgb "red" lw 2
set linestyle 3 lt 6 lc rgb "orange" lw 2

set border 3
set tics nomirror
set grid
set key r b
set xrange [0:180]
set xtics 20
set xlabel "Rotation angle (degrees)"

s4fForce = "postProcessing/0/solidForcessphere-displacement.dat"
s4fTorque = "postProcessing/0/solidTorquesphere-displacementsphereTorque.dat"

# Vertical force
set output "force.png"
set ylabel "Vertical force"
set yrange [0:1.5]
set ytics 0.25

plot \
    s4fForce u (20*($1 - 4)):(-$3) w lp ls 1 pi 5 title "solids4foam", \
    "deLorenzisForce.dat" u 1:2 w l ls 2 title "Sauer and De Lorenzis"

# Twisting torque
set output "torque.png"
set ylabel "Twisting torque"
set yrange [0:0.5]
set ytics 0.1

plot \
    s4fTorque u (20*($1 - 4)):2 w lp ls 1 pi 5 title "solids4foam", \
    "deLorenzisMoment.dat" u 1:2 w l ls 2 title "Sauer and De Lorenzis", \
    "febioMoment.dat" u 1:($2/10) w l ls 3 title "Zimmerman and Ateshian (FEBio)"
