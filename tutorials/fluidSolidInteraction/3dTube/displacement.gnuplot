set terminal pdfcairo enhanced color solid

# Axial displacement

set output "axialDisplacement.pdf"
set title "3dTube: Axial Displacement at Point A"
set xlabel "Time, t [s]"
set ylabel "Axial Displacement [mm]"
set grid
set key top right

set yrange [-0.1:0.1]

plot \
    "postProcessing/0/solidPointDisplacement_displacement.dat" u 1:(1000*$4)\
    w l lc "red" lw 0.5 lt 1 t "Axial displacement"

# Radial displacement

set output "radialDisplacement.pdf"
set title "3dTube: Radial Displacement at Point A"
set ylabel "Radial Displacement [mm]"
set yrange [-0.1:0.15]

plot \
    "postProcessing/0/solidPointDisplacement_displacement.dat" u 1:(1000*$3)\
    w l lc "red" lw 0.5 lt 1 t "Radial displacement"
