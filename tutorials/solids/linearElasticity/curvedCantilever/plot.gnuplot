#!/usr/bin/gnuplot

set term pngcairo dashed font "Arial,20" size 1920,1080
set output 'sigmaAtTheta45deg.png'

set linestyle  1 lt 6 lc rgb "black" lw 1.5 ps 1 dt 1
set linestyle  2 lt 6 lc rgb "black" lw 1.5 ps 1 dt 2
set linestyle  3 lt 6 lc rgb "black" lw 1.5 ps 1 dt 3
set linestyle  4 lt 6 lc rgb "red" lw 1.5 ps 1 dt 1
set linestyle  5 lt 6 lc rgb "red" lw 1.5 ps 1 dt 2
set linestyle  6 lt 6 lc rgb "red" lw 1.5 ps 1 dt 3

set grid
set title "Stress distribution along line at {/Symbol q} = 45deg" font "Arial Bold,20"
set xlabel "x [mm]" font "Arial Bold,20"
set ylabel "Stress [Pa]" font "Arial Bold,20"
set key r b

plot \
"postProcessing/sets/1/line_sigma_analyticalStress.xy" u 1:4 w lp ls 1 title"{/Symbol s}_{xx}",\
"postProcessing/sets/1/line_sigma_analyticalStress.xy" u 1:7 w lp ls 2 title"{/Symbol s}_{yy}",\
"postProcessing/sets/1/line_sigma_analyticalStress.xy" u 1:5 w lp ls 3 title"{/Symbol s}_{xy}",\
"postProcessing/sets/1/line_sigma_analyticalStress.xy" u 1:10 w lp ls 4 title"Analytical {/Symbol s}_{xx}",\
"postProcessing/sets/1/line_sigma_analyticalStress.xy" u 1:13 w lp ls 5 title"Analytical {/Symbol s}_{yy}",\
"postProcessing/sets/1/line_sigma_analyticalStress.xy" u 1:11 w lp ls 6 title"Analytical {/Symbol s}_{xy}"
