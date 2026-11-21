#
# Analytical symbols (lines) and measured eigenvalues (points) of the raw
# pressure-stabilisation operators along the cut theta_y = 0, at 16 cells per
# axis, each scaled by its analytical value at theta_x = pi
#

reset
set term pngcairo dashed size 1024,768 font "Arial,16"
dir = "postProcessing/stabilisationFourierCheck"
set output dir . "/symbols_N16.png"

set datafile separator ","
set grid
set key left top
set yrange [0:1.05]
set xlabel "theta_x / pi"
set ylabel "lambda / lambda_{exact}(theta_x = pi)"

models = "laplacian gen1 gen2 diffStencil"

# Analytical cut 0 (theta_y = 0) from symbols_*.csv, scaled by its last value
cut(m) = sprintf("< awk -F, '$1 == \"%s\" && $2 == 0 {print $3\",\"$4}' %s/symbols_16x16x1.csv", m, dir)
ref(m) = sprintf("awk -F, '$1 == \"%s\" && $2 == 0 {v = $4} END {print v}' %s/symbols_16x16x1.csv", m, dir)
# Measured raw eigenvalues on the same cut from results_*.csv
# (columns: normalise, model, thetaX/pi, thetaY/pi, thetaZ/pi, lambda_num)
pts(m) = sprintf("< awk -F, '$1 == 0 && $3 == \"%s\" && $21 == 0 && $22 == 0 {print $20\",\"$24}' %s/results_16x16x1.csv", m, dir)

plot for [i=1:words(models)] cut(word(models, i)) \
         using 1:($2/real(system(ref(word(models, i))))) \
         with lines lw 2 lc i title word(models, i), \
     for [i=1:words(models)] pts(word(models, i)) \
         using 1:($2/real(system(ref(word(models, i))))) \
         with points pt 7 ps 1.5 lc i notitle
