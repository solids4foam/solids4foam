# Drag and lift coefficients of the translating immersed cylinder, compared
# with those of the static cylinder at the same positions
# (verificationData/CdStatic.dat)
#
# Cd = 2 Fx/(rho Umean^2 D Lz), with rho = 1 kg/m^3, the mean velocity of the
# flow relative to the cylinder Umean = 0.2 m/s, D = 0.1 m and Lz = 0.1 m,
# i.e. Cd = 5000 Fx. Columns 2-4 of the force file are the force from the
# momentum exchange, and columns 11-13 that from the surface traction.
set terminal pdfcairo enhanced color font "Helvetica,12" size 6,6
set output "forceCoeffs.pdf"

forceFile = "postProcessing/immersedBoundary/0/cylinder.dat"
refFile = "verificationData/CdStatic.dat"

set multiplot layout 2,1

set xlabel "Time (s)"
set ylabel "C_d"
set xrange [1:10]
set yrange [4.8:6]
set grid
plot \
    forceFile u 1:(5000*$2) w l lw 2 lc rgb "black" \
        t "immersedBoundaryForce (momentum exchange)", \
    forceFile u 1:(5000*$11) w l lw 1 dt 2 lc rgb "blue" \
        t "immersedBoundaryForce (surface traction)", \
    refFile u 1:3 w lp pt 7 ps 0.5 lc rgb "red" \
        t "Static cylinder at the same position"

set ylabel "C_l"
set yrange [-0.05:0.05]
plot \
    forceFile u 1:(5000*$3) w l lw 2 lc rgb "black" \
        t "immersedBoundaryForce (momentum exchange)", \
    refFile u 1:4 w lp pt 7 ps 0.5 lc rgb "red" \
        t "Static cylinder at the same position"

unset multiplot
