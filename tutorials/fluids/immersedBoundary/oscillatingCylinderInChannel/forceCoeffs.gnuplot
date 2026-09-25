# Drag and lift coefficients of the immersed oscillating cylinder, compared
# with Wan and Turek (2006) and a moving body-fitted mesh solution
#
# C = 2 F/(rho Uref^2 D Lz), with rho = 1 kg/m^3, the maximum cylinder
# velocity Uref = 2 pi A/T = 0.3927 m/s, D = 0.1 m and Lz = 0.1 m, i.e.
# C = 1296.9 F. The hydrodynamic force is the force of the forcing (columns
# 2-4) plus the inertia of the fluid inside the cylinder (columns 8-10).
set terminal pdfcairo enhanced color font "Helvetica,12" size 6,6
set output "forceCoeffs.pdf"

forceFile = "postProcessing/immersedBoundary/0/cylinder.dat"
scale = 1296.9

set multiplot layout 2,1

set xlabel "Time (s)"
set ylabel "C_d"
set yrange [-6:6]
set grid
plot \
    forceFile u 1:(scale*($2 + $8)) w l lw 2 lc rgb "black" \
        t "immersedBoundaryForce", \
    "verificationData/Cd.dat" u 1:2 w p pt 7 ps 0.3 lc rgb "red" \
        t "Wan and Turek (2006)", \
    "verificationData/CdBodyFitted.dat" u 1:2 w l lw 1 dt 2 lc rgb "blue" \
        t "Body-fitted mesh"

set ylabel "C_l"
set yrange [-0.1:0.1]
plot \
    forceFile u 1:(scale*($3 + $9)) w l lw 2 lc rgb "black" \
        t "immersedBoundaryForce", \
    "verificationData/Cl.dat" u 1:2 w p pt 7 ps 0.3 lc rgb "red" \
        t "Wan and Turek (2006)", \
    "verificationData/ClBodyFitted.dat" u 1:2 w l lw 1 dt 2 lc rgb "blue" \
        t "Body-fitted mesh"

unset multiplot
