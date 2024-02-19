reset
set term pngcairo dashed size 1024,768 font "Helvetica,18"

set style line 1 linecolor rgb 'black' linetype 6 linewidth 2 ps 1

set grid
set xlabel "Normalised Strain (dimensionless)"
set ylabel "Normalised Stress (in MPa)"
set key right top
set xrange [0:0.005]
set yrange [0:360]
set xtics 0.001
set ytics 50
set key top left

set output 'stressVsStrain.png'

plot "postProcessing/0/solidForcesup.dat" u ($1*2*1.143/152.4):(2*$3/((25.4 - 2*4.06)*2.37)) w l ls 1 title "solids4foam"
