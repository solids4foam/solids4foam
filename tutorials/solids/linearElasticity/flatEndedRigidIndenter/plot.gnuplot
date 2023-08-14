#!/usr/bin/gnuplot

reset
set term pngcairo dashed font "Arial,20" size 1520,1080
set output "pressureDistribution.png"

set linestyle  1 lt 6 lc rgb "red" lw 1.5 ps 2
set linestyle  2 lt 6 lc rgb "black" lw 1.5 ps 2

set datafile separator " "
set xrange [-0.015:0.015]
set yrange [0:1500]
set key center top
set grid
set xlabel "x [m]" font "Arial,20"
set ylabel "p_n [MPa]" font "Arial,20"

path = "./postProcessing/surfaces/1/sigma_surface.raw"

a = 10    # In mm to get MPa
F = 10000 # In N
pn(x) = (x > -0.00995 && x < 0.00995) ? F/(pi*(a**2-(x*1e3)**2)**(0.5)) : 1/0
set samples 1000000

plot [-0.015:0.015] pn(x) w l ls 1 title"Analytical {p}_{n}(x)",\
     path u 1:(abs($7)*1e-6) w lp ls 2 notitle smooth unique,\
     path u 1:(abs($7)*1e-6) w p ls 2 notitle
