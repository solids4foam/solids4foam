# Drag and lift coefficients of the bending immersed beam, compared with those
# of the body-fitted solution with a deforming mesh
# (verificationData/CdBodyFitted.dat)
#
# Cd = 2 Fx/(rho Umean^2 H Lz), with rho = 1 kg/m^3, the mean inlet velocity
# Umean = 1 m/s, the beam height H = 2.0295 m and Lz = 0.1 m, i.e.
# Cd = 9.855 Fx. Columns 2-4 of the force file are the force from the surface
# traction, and columns 11-13 that from the momentum exchange.
set terminal pdfcairo enhanced color font "Helvetica,12" size 6,6
set output "forceCoeffs.pdf"

forceFile = "postProcessing/immersedBoundary/0/beam.dat"
refFile = "verificationData/CdBodyFitted.dat"

set multiplot layout 2,1

set xlabel "Time (s)"
set ylabel "C_d"
set xrange [0:8]
set yrange [-80:120]
set grid
set key top left horizontal
plot \
    forceFile u 1:(9.855*$2) w l lw 2 lc rgb "black" \
        t "Surface traction", \
    forceFile u 1:(9.855*$11) w l lw 1 dt 2 lc rgb "blue" \
        t "Momentum exchange", \
    refFile u 1:2 w l lw 2 lc rgb "red" \
        t "Body-fitted"

set ylabel "C_l"
set yrange [-20:25]
plot \
    forceFile u 1:(9.855*$3) w l lw 2 lc rgb "black" \
        t "Surface traction", \
    refFile u 1:3 w l lw 2 lc rgb "red" \
        t "Body-fitted"

unset multiplot
