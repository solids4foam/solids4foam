# Deformed interface at t = 1 s and in the steady state, compared with Liu
# (arXiv:1401.0082). Called by blob_in_treacle_verification.py with the
# variables 'interfaces', 'reference' and 'output'.

set terminal pngcairo enhanced size 1400,560 font ",12"
set output output
set datafile separator ","
set key bottom center
set grid
set size ratio -1
set xlabel "x (m)"
set ylabel "y (m)"
set xrange [0.95:2.1]
set yrange [-0.52:0.05]

set multiplot layout 1,2

set title "t = 1 s"
plot \
    reference."/LiuInterface_t1.csv" using 1:2 with points pt 7 ps 0.3 lc rgb "black" title "Liu, dt = 0.025 s", \
    interfaces."/mesh_x1_t1.csv" using 1:2 with linespoints pt 6 ps 0.6 lc rgb "#1f77b4" title "solids4foam, tutorial mesh", \
    for [f in "2 4"] interfaces."/mesh_x".f."_t1.csv" using 1:2 with lines lw 2 title "solids4foam, mesh x".f

set title "Steady state (Liu: t = 10 s; solids4foam: t = 3 s)"
plot \
    reference."/LiuInterface_t10.csv" using 1:2 with points pt 7 ps 0.8 lc rgb "black" title "Liu, domain decomposition", \
    interfaces."/mesh_x1_steady.csv" using 1:2 with linespoints pt 6 ps 0.6 lc rgb "#1f77b4" title "solids4foam, tutorial mesh", \
    for [f in "2 4"] interfaces."/mesh_x".f."_steady.csv" using 1:2 with lines lw 2 title "solids4foam, mesh x".f

unset multiplot
