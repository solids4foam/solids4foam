reset
set term pngcairo dashed size 1024,768 font "Arial,18"

set style line 1 linecolor rgb 'black' linetype 6 linewidth 2 ps 1.8
set style line 2 linecolor rgb 'red' linetype 6 linewidth 2 ps 1.8

set grid
set xlabel "Normalised displacement"
set ylabel "Reaction-x"
set key left top
set xrange [0:1]
set yrange [-30:10]
set xtics 0.2
set output 'reaction-x.png'

plot "postProcessing/0/solidForcesfixed.dat" u ($1/31.5):($2) w lp ls 1 title"solids4foam",\
     "netoEtAl.dat" u ($1/31.5):2 w p ls 2 title"Neto et al."
     
set output 'reaction-y.png'
set yrange [0:50]
set ylabel "Reaction-y"

plot "postProcessing/0/solidForcesfixed.dat" u ($1/31.5):($3) w lp ls 1 title"solids4foam",\
     "netoEtAl.dat" u ($3/31.5):4 w p ls 2 title"Neto et al."

