reset
set term pngcairo dashed size 1024,768 font "Helvetica,18"

set style line 1 linecolor rgb 'black' linetype 1 linewidth 2
set style line 2 linecolor rgb 'red' linetype 2 linewidth 2 dashtype 2

set grid
set xlabel "Displacement (in mm)"
set ylabel "Force (in kN)"
set key right bottom
set xrange [0:1]
set yrange [0:40]

set output 'forceDisplacement.png'

# The force on the 5 degree wedge is scaled by 360/5 to the whole bar
set datafile separator whitespace
plot "postProcessing/0/solidForcesDisplacementsup.dat" \
    u ($3*1000):($6*72e-3) w l ls 1 title "solids4foam", \
    "< tr ',' ' ' < reference/abaqusForceDisplacement.csv" \
    u 1:2 w l ls 2 title "Abaqus"
