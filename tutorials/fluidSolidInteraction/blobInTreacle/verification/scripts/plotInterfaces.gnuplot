# Deformed interface at t = 1 s and in the steady state, compared with Liu
# (arXiv:1401.0082). Called by blob_in_treacle_verification.py with the
# variables 'interfaces', 'reference', 'output' and 'levels', the list of
# completed mesh refinement factors, e.g. "0.5 1 2 4".

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

# Level 1 is the tutorial mesh
label(f) = (f eq "1") ? "solids4foam, tutorial mesh" : "solids4foam, mesh x".f

set multiplot layout 1,2

set title "t = 1 s"
plot \
    reference."/LiuInterface_t1.csv" using 1:2 with points pt 7 ps 0.3 lc rgb "black" title "Liu, dt = 0.025 s", \
    for [f in levels] interfaces."/mesh_x".f."_t1.csv" using 1:2 with lines lw 2 title label(f)

set title "Steady state (Liu: t = 10 s; solids4foam: t = 3 s)"
plot \
    reference."/LiuInterface_t10.csv" using 1:2 with points pt 7 ps 0.8 lc rgb "black" title "Liu, domain decomposition", \
    for [f in levels] interfaces."/mesh_x".f."_steady.csv" using 1:2 with lines lw 2 title label(f)

unset multiplot
