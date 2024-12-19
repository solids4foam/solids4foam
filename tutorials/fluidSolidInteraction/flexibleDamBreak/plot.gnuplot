set terminal pdfcairo enhanced color solid

set output "displacement.pdf"

set style line 1 linecolor rgb 'black' linetype 6 linewidth 2 ps 1.8
set style line 2 linecolor rgb 'red' linetype 6 linewidth 2 ps 1.8
set style line 3 linecolor rgb 'blue' linetype 6 linewidth 2 ps 1.8
set style line 4 linecolor rgb 'green' linetype 6 linewidth 2 ps 1.8
set key t r
set key box
set xrange [0:1]
set yrange [-0.03:0.055]
set xtic nomirror
set ytic nomirror
set mxtic 4
set mytic 4
plot "postProcessing/0/solidPointDisplacement_disp.dat" u ($1):($2) w l ls 1 title "solids4foam" \
        ,"referenceSolution/refSol.tsv" u ($1):($2) w l ls 2 title "Meduri et al. [2]" \
        ,"referenceSolution/refSol.tsv" u ($1):($3) w l ls 3 title "Ryzhakov et al. [3]"
