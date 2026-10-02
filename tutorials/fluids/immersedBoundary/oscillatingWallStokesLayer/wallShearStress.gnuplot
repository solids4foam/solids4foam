# Wall shear stress on the immersed oscillating slab, from the force on the
# slab divided by the wall area in the mesh (0.05 m x 0.1 m), compared with
# the periodic solution of the Stokes second problem:
#     tau = -nu k U0 (sin(omega t) + cos(omega t)), k = sqrt(omega/(2 nu))
# with U0 = 0.1 m/s, omega = 2 pi rad/s and nu = 0.00785398 m^2/s (rho = 1)
set terminal pdfcairo enhanced color font "Helvetica,12" size 6,4
set output "wallShearStress.pdf"

forceFile = "postProcessing/immersedBoundary/0/slab.dat"
A = 0.005
nu = 0.00785398163
U0 = 0.1
omega = 2*pi
k = sqrt(omega/(2*nu))
tau(t) = -nu*k*U0*(sin(omega*t) + cos(omega*t))

set xlabel "Time (s)"
set ylabel "Wall shear stress (Pa)"
set xrange [2:4]
set grid
set samples 400
plot \
    forceFile u 1:(($2 + $8)/A) w l lw 2 lc rgb "black" \
        t "immersedBoundaryForce (momentum exchange)", \
    forceFile u 1:($11/A) w l lw 1 dt 2 lc rgb "blue" \
        t "immersedBoundaryForce (surface traction)", \
    tau(x) w l lw 1 lc rgb "red" t "Periodic solution"
