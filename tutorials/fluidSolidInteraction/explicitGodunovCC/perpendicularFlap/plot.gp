# Gnuplot script to plot Time vs Dx
set terminal pdfcairo enhanced color solid

set output "displacement.pdf"

# Set labels
set xlabel 'Time (s)'
set ylabel 'Dx (m)'
set title 'Time vs Dx Displacement'

set style line 1 linecolor rgb 'black' linetype 6 linewidth 2 ps 1.8
set style line 2 linecolor rgb 'red' linetype 6 linewidth 2 ps 1.8
set style line 3 linecolor rgb 'blue' linetype 6 linewidth 2 ps 1.8
set style line 4 linecolor rgb 'green' linetype 6 linewidth 2 ps 1.8
set key t r
set key box
set xrange [0:5]
set yrange [-0.03:0.25]
set xtic nomirror
set ytic nomirror
set mxtic 4
set mytic 4

# Plotting the data
plot 'postProcessing/0/solidPointDisplacement_displacement.dat' u ($1):($2) w l ls 1  title 'Dx Displacement'\
    ,"referenceSolution/refSol.dat" u ($1):($2) w l ls 2 title "preCICE. [2]" \

