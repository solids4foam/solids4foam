# Displacements of the benchmark points against the benchmark participants'
# solutions, for the ArosticaEtAl2025 tutorials. Run from the case directory:
#
#     gnuplot -e "refDir='verification/reference/step1'; points='p0 p1'; \
#         output='verification/step1.pdf'" ../plotVsReference.gnuplot
#
# refDir holds one CSV per participant solution (time, then ux, uy, uz of each
# point, in m), from the cardiac-benchmark-toolkit (MIT licence). The
# solids4foam histories are postProcessing/0/solidPointDisplacement_<p>.dat
# (time, ux, uy, uz, |u|, in m). One panel per point and component, in mm.
#
# Created by Claude Opus 5.5, supervised by Chanikya Valeti, UCD.

set terminal pdfcairo enhanced color size 5*words(points),9 font ",9"
set output output

refFiles = system("ls ".refDir."/*.csv")
nP = words(points)
comps = "x y z"

set multiplot layout 3, nP title "Aróstica et al. (2025): solids4foam vs participants (".refDir.")"
set grid
set key top left box opaque

do for [c = 1:3] {
    do for [i = 1:nP] {
        p = word(points, i)
        set title p." u_".word(comps, c)
        set xlabel (c == 3 ? "t [s]" : "")
        set ylabel "displacement [mm]"
        refCol = 1 + 3*(i - 1) + c
        s4fFile = "postProcessing/0/solidPointDisplacement_".p.".dat"
        plot for [j = 1:words(refFiles)] \
                "< tr ',' ' ' < ".word(refFiles, j) skip 1 \
                using 1:(1e3*column(refCol)) with lines lc rgb "#b0b0b0" \
                lw 1 title (j == 1 ? "participants" : ""), \
             s4fFile using 1:(1e3*column(c + 1)) with lines lc rgb "#d62728" \
                lw 2 title "solids4foam"
    }
}

unset multiplot
