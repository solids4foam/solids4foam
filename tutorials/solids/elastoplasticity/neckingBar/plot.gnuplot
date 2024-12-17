reset
set term pngcairo dashed size 1024,768 font "Helvetica,18"

set style line 1 linecolor rgb 'black' linetype 6 linewidth 2 ps 1

set grid
set xlabel "Axial elongation (in mm)"
set ylabel "Force (in kN)"
set key right top
set xrange [0:7]
set yrange [0:80]
set xtics 1
set ytics 10

set output 'forceElongation.png'

plot "postProcessing/0/solidForcesDisplacementsloading.dat" u ($2*1000):($5*360*1e-3) w l ls 1 title"solids4foam"
