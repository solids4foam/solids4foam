# Plot the radial and axial wall displacement at point A and the axis
# pressure at z = 2.5 cm for a set of verification runs, with any published
# point-A histories overlaid.
#
# Invoked by the verification driver as:
#   gnuplot -e "files='a.csv b.csv'; titles='a b'; references='r.csv';
#                referenceTitles='r'; axialReferences='r.csv';
#                axialReferenceTitles='r'; term='png'; output='...'"
#                plotPointAHistory.gnuplot

if (term eq "png") {
    set terminal pngcairo enhanced size 900,1100
} else {
    set terminal pdfcairo enhanced size 6in,7.5in
}
set output output

set datafile separator ","
set datafile missing "nan"
set multiplot layout 3,1
set grid
set key top right
set xlabel "Time (ms)"
set xrange [0:20]

set ylabel "Radial displacement at A (mm)"
plot for [i=1:words(files)] word(files, i) using ($1*1e3):($2*1e3) \
         with lines lw 1.5 title word(titles, i) noenhanced, \
     for [i=1:words(references)] word(references, i) \
         using ($1*1e3):($2*1e3) with lines lw 1.5 dt 2 \
         title word(referenceTitles, i) noenhanced

set ylabel "Axial displacement at A (mm)"
plot for [i=1:words(files)] word(files, i) using ($1*1e3):($3*1e3) \
         with lines lw 1.5 title word(titles, i) noenhanced, \
     for [i=1:words(axialReferences)] word(axialReferences, i) \
         using ($1*1e3):($3*1e3) with lines lw 1.5 dt 2 \
         title word(axialReferenceTitles, i) noenhanced

set ylabel "Axis pressure at z = 2.5 cm (Pa)"
plot for [i=1:words(files)] word(files, i) using ($1*1e3):4 \
         with lines lw 1.5 title word(titles, i) noenhanced

unset multiplot
