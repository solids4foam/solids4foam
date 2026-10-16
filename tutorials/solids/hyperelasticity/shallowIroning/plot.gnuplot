reset
set term pngcairo dashed size 1024,768 font "Arial,16"

set style line 1 linecolor rgb 'black' linetype 1 linewidth 3
set style line 2 linecolor rgb 'red' linetype 1 linewidth 1
set style line 3 linecolor rgb 'blue' linetype 1 linewidth 1
set style line 4 linecolor rgb 'dark-green' linetype 1 linewidth 1
set style line 5 linecolor rgb 'dark-cyan' linetype 1 linewidth 1
set style line 6 linecolor rgb 'orange' linetype 1 linewidth 1

# The solids4foam time t in (1, 6] s is mapped to the normalised time
# 1 + (t - 1)/5 in (1, 2] used by the reference data. The reference forces are
# multiplied by 100 to compare with the solids4foam force (N) on the 1 mm thick
# domain.
tn(t) = (t > 1) ? 1 + (t - 1)/5.0 : t

set grid
set xlabel "Normalised time"
set ylabel "Reaction force (N)"
set key left top
set xrange [0:2]
set yrange [0:450]
set output 'reactionForces.png'

plot "postProcessing/0/solidForcesdisplacement.dat" u (tn($1)):(-$3) w l ls 1 title "solids4foam (vertical)", \
     "postProcessing/0/solidForcesdisplacement.dat" u (tn($1)):($2) w l ls 1 dt 2 title "solids4foam (horizontal)", \
     "reference/asterVertical.dat" u 1:(100*$2) w l ls 2 title "Code\\_Aster", \
     "reference/asterHorizontal.dat" u 1:(100*$2) w l ls 2 notitle, \
     "reference/hartmannOliverVertical.dat" u 1:(100*$2) w l ls 3 title "Hartmann et al.", \
     "reference/hartmannOliverHorizontal.dat" u 1:(100*$2) w l ls 3 notitle, \
     "reference/fischerWriggersVertical.dat" u 1:(100*$2) w l ls 4 title "Fischer and Wriggers", \
     "reference/fischerWriggersHorizontal.dat" u 1:(100*$2) w l ls 4 notitle, \
     "reference/pouliosRenardVertical.dat" u 1:(100*$2) w l ls 5 title "Poulios and Renard", \
     "reference/pouliosRenardHorizontal.dat" u 1:(100*$2) w l ls 5 notitle, \
     "reference/yastrebov.dat" u 1:(-100*$3) w l ls 6 title "Yastrebov", \
     "reference/yastrebov.dat" u 1:(100*$2) w l ls 6 notitle
