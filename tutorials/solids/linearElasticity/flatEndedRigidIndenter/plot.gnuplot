#!/usr/bin/gnuplot

reset
set term pngcairo dashed font "Helvetica,20" size 1520,1080
set output "pressureDistribution.png"

set linestyle  1 lt 6 lc rgb "red" lw 1.5 ps 2
set linestyle  2 lt 6 lc rgb "black" lw 1.5 ps 2

# Analytical solution
a = 0.01   # Contact half-width (in m)
F = 10000 # Applied force on indenter (in N)
d = -0.09673053732e-3 # Indenter displacement (in m)
nu = 0.3 # Poissons ratio
E = 200000e6 # Young modulus (in Pa)

set samples 1000000
pn(x) = (x > -0.00995 && x < 0.00995) ? F/(pi*sqrt((a*1e3)**2-(x*1e3)**2)) : 1/0
u(x) = (x > 0.01 ) ? 1e3*(d+1000*((2.0*(1.0-nu**2)*F)/(pi*E))*log((x/a)+sqrt((x/a)**2-1.0))) : 1/0
d(x) = (x < 0.01 ) ? 1e3*d : 1/0

set datafile separator " "
set yrange [0:1500]
set key center top
set grid
set xlabel "x [m]"
set ylabel "p_n [MPa]"

pathSigma = "postProcessing/surfaces/1/sigma_surface.raw"

plot [-0.015:0.015] pn(x) w l ls 1 title"Analytical",\
     pathSigma u 1:(abs($7)*1e-6) w lp ls 2 notitle smooth unique,\
     pathSigma u 1:(abs($7)*1e-6) w p ls 2 title"solids4foam"

set output "displacement.png"
set xrange [0:0.05]
set yrange [-0.1:-0.02]
set key b r
set ylabel "u_z [mm]"

pathD = "postProcessing/surfaces/1/D_surface.raw"

plot u(x) w l ls 1 title"Analytical", d(x) w l ls 1 notitle,\
     pathD u 1:($5*1e3) w p ls 2 title"solids4foam"
