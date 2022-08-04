set xlabel "Time (in s)"
set ylabel "Energy (in J)"

set style line 81 lt 0 lw 0.5
set style line 81 lt rgb "#808080"
set grid back linestyle 81
set border 3 back linestyle 80

set xtics nomirror
set ytics nomirror

set title "Energies"

plot \
    "energies.dat" u 1:2 w l lw 2 t "External work",\
    "" u 1:3 w l lw 2 t "Internal/strain (I)", \
    "" u 1:4 w l lw 2 t "Kinetic (K)", \
    "" u 1:5 w l lw 2 t "Viscous (V)", \
    "" u 1:($3 + $4 + $5) w l lw 2 t "Sum(I+K+V)"

pause 1
reread