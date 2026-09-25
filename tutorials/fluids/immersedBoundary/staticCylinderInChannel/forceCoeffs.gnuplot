# Drag and lift coefficients of the immersed cylinder, compared with the
# bounds of Schafer and Turek (1996)
#
# Cd = 2 Fx/(rho Umean^2 D Lz), with rho = 1 kg/m^3, Umean = 0.2 m/s,
# D = 0.1 m and Lz = 0.1 m, i.e. Cd = 5000 Fx
set terminal pdfcairo enhanced color font "Helvetica,12" size 6,6
set output "forceCoeffs.pdf"

forceFile = "postProcessing/immersedBoundary/0/cylinder.dat"

set multiplot layout 2,1

set xlabel "Time (s)"
set ylabel "C_d"
set yrange [5:7]
set grid
plot \
    forceFile u 1:(5000*$2) w l lw 2 lc rgb "black" t "immersedBoundaryForce", \
    5.57 w l dt 2 lc rgb "red" t "Schafer and Turek (1996) bounds", \
    5.59 w l dt 2 lc rgb "red" notitle

set ylabel "C_l"
set yrange [-0.05:0.05]
plot \
    forceFile u 1:(5000*$3) w l lw 2 lc rgb "black" t "immersedBoundaryForce", \
    0.0104 w l dt 2 lc rgb "red" t "Schafer and Turek (1996) bounds", \
    0.0110 w l dt 2 lc rgb "red" notitle

unset multiplot
