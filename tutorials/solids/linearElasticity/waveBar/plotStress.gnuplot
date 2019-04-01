set xlabel "Time (in s)"
set ylabel "Stress (in MPa)"

set style line 81 lt 0 lw 0.5
set style line 81 lt rgb "#808080"
set grid back linestyle 81
set border 3 back linestyle 80

set xtics nomirror
set ytics nomirror

plot \
    "< sed s/[\\(\\)]//g probe/0/sigma" using 1:(1e-6*$2) lw 2 lt 1 w l \
    t "Axial stress"