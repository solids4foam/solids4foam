# Compression force versus displacement of the upper rigid plane
#
# Usage: gnuplot -e "contact='frictionless'" plot.gnuplot
#    or: gnuplot -e "contact='friction'" plot.gnuplot
# The contact type defaults to frictionless

if (!exists("contact")) contact = "frictionless"

set term pngcairo dashed size 1024,768 font "Arial,18"

set style line 1 linecolor rgb "black" linetype 1 linewidth 2 pointtype 7 ps 0.8
set style line 2 linecolor rgb "red" linetype 1 linewidth 2
set style line 3 linecolor rgb "orange" linetype 1 linewidth 2
set style line 4 linecolor rgb "forest-green" linetype 1 linewidth 2
set style line 5 linecolor rgb "purple" linetype 1 linewidth 2 pointtype 6 ps 1.2

set border 3
set tics nomirror
set grid
set xlabel "Displacement of the upper plane (mm)"
set ylabel "Compression force (N)"
set xrange [0:10]
set yrange [0:14]

s4f = "postProcessing/0/solidForcesR_top.dat"
ref = "reference/".contact."/"

set output "force-displacement.png"

# The solids4foam time goes from 0 to 1 while the plate moves 10 mm; the
# reference data unit conversions are described in the reference data headers
if (contact eq "friction") {
    set title "Frictional contact ({/Symbol m} = 0.5)"
    set key left top
    plot ref."Abaqus.dat" u 1:2 w l ls 4 title "Abaqus [2]", \
         ref."FEBio.dat" u 1:2 w l ls 3 title "FEBio [2]", \
         ref."Areias.dat" u ($1*10):(-$2*100) w l ls 2 title "Areias et al. [4]", \
         s4f u ($1*10):4 w lp ls 1 title "solids4foam"
} else {
    set title "Frictionless contact"
    set key left top
    plot ref."Abaqus.dat" u 1:2 w l ls 4 title "Abaqus [2]", \
         ref."FEBio.dat" u ($1*10):(-$2) w l ls 3 title "FEBio [2]", \
         ref."Areias.dat" u ($1*10):(-$2) w l ls 2 title "Areias et al. [4]", \
         ref."PusoLaursen.dat" u ($1*10):(-$2) w p ls 5 title "Puso and Laursen [3]", \
         s4f u ($1*10):4 w lp ls 1 title "solids4foam"
}
