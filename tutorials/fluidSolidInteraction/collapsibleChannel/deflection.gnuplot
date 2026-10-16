set terminal pdfcairo enhanced color solid
set output "deflection.pdf"

set xlabel "Time (s)"
set ylabel "Vertical wall displacement (m)"
set grid
set key bottom right

plot \
    "postProcessing/0/solidPointDisplacement_wallQuarter.dat" \
        using 1:3 with lines lw 2 title "x = 7.5 m", \
    "postProcessing/0/solidPointDisplacement_wallMid.dat" \
        using 1:3 with lines lw 2 title "x = 10 m", \
    "postProcessing/0/solidPointDisplacement_wallThreeQuarter.dat" \
        using 1:3 with lines lw 2 title "x = 12.5 m"
