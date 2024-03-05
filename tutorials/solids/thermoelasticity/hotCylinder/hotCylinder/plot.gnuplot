#!/usr/bin/gnuplot

set term pngcairo dashed font "Times,26" size 1520,1080

set grid
set xrange [0.5:0.7]
set xtics 0.05
set xlabel "r [m]" font "Times,20"
set key r t
set datafile separator " "

set style line 1 linecolor rgb 'red' linetype 6 linewidth 2 ps 2
set style line 2 linecolor rgb 'black' linetype 6 linewidth 2 ps 2

path2="postProcessing/sets/1000/line_sigma\:Transformed.xy"
path1="postProcessing/sets/1000/line_analyticalRadialStress_analyticalHoopStress_T_analyticalT.xy"

set output 'sigmaR.png'
set yrange [0:-15]
set ytics 2.5
set ylabel "{/Symbol s}_{r} [MPa]"
plot path2 using 1:($2*1e-6) w lp ls 2 title"solids4foam",\
     path1 using 1:($2*1e-6) w lp ls 1 title"Analytical"

set output 'sigmaTheta.png'
set yrange [-200:250]
set ytics 50
set ylabel "{/Symbol s}_{{/Symbol q}} [MPa]"
plot path2 using 1:($5*1e-6) w lp ls 2 title"solids4foam",\
     path1 using 1:($3*1e-6) w lp ls 1 title"Analytical"
     
set output 'T.png'
set yrange [0:100]
set ytics 10
set ylabel "T [K]"
plot path1 using 1:4 w lp ls 2 title"solids4foam",\
     path1 using 1:5 w lp ls 1 title"Analytical"
 

