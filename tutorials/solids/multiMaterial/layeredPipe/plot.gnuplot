#!/usr/bin/gnuplot

set term pngcairo dashed font "Times,20" size 1920,1080

set grid
set xrange [50:100]
set xtics 5
set xlabel "r [mm]"
set ylabel "Stress [in MPa]" font "Times Bold,20"
set key outside center top horizontal

set output 'sigmaXX_lineXX.png'
set yrange [0:0.1]
set ytics 0.01
plot "postProcessing/sample/10/lineXX_sigma.xy" using ($1*1000):(-$2*1e-6) w l

set output 'sigmaXX_lineYY.png'
set yrange [-0.05:0.25]
set ytics 0.05
plot "postProcessing/sample/10/lineYY_sigma.xy" using ($1*1000):($2*1e-6) w l
