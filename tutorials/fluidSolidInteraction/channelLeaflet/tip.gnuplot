set terminal pdfcairo enhanced color solid
set output "tip.pdf"

set xlabel "Time (s)"
set ylabel "Leaflet tip displacement (m)"
set grid
set key bottom left

plot \
    "postProcessing/0/solidPointDisplacement_tip.dat" \
        using 1:2 with lines lw 2 title "x", \
    "postProcessing/0/solidPointDisplacement_tip.dat" \
        using 1:3 with lines lw 2 title "y"
