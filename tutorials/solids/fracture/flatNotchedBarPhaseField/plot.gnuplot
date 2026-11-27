reset
set term pngcairo dashed size 1024,768 font "Helvetica,18"

set style line 1 linecolor rgb 'black' linetype 1 linewidth 2
set style line 2 linecolor rgb 'red' linetype 2 linewidth 2 dashtype 2
set style line 3 linecolor rgb 'blue' pointtype 7 pointsize 1.2

set grid
set xlabel "Normalised displacement (-)"
set ylabel "Normalised stress (in MPa)"
set key right bottom
set xrange [0:0.005]
set yrange [0:400]

set output 'stressDisplacement.png'

# The force on the quarter plate over the area of its notched section,
# 8.64 mm by 2.37 mm, against the displacement over the gauge half-length,
# 76.2 mm
set datafile separator whitespace
plot "postProcessing/0/solidForcesDisplacementsup.dat" \
    u ($3/76.2e-3):($6/20.4768) w l ls 1 title "solids4foam", \
    "< tail -n +2 reference/bordenNormalisedStressStrain.csv | tr ',' ' '" \
    u 1:2 w l ls 2 title "Borden et al. (2016)", \
    "< tail -n +2 reference/eldahshanNormalisedStressStrain.csv | tr ',' ' '" \
    u 1:2 w p ls 3 title "Eldahshan et al. (2021)"
