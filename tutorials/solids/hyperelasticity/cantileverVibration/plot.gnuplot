# Single selected run, or all three time schemes from ./Allrun petscSnes all
set term pngcairo dashed size 1200,850 font "Arial,18"
set output "tipDisplacement.png"
set grid
set xlabel "Time (s)"
set ylabel "Tip displacement magnitude (m)"
set key right top
set xrange [0:0.65]
if (!exists("comparison")) comparison = 0
if (!exists("scheme")) scheme = "bdf2"

if (comparison) {
    plot \
        "reference/abaqusC3D8.dat" u 1:2 w l lw 3 lc "black" t "Abaqus (C3D8)", \
        "timeSchemeRuns/bdf2/postProcessing/0/solidPointDisplacement_pointDisp.dat" \
            u 1:5 w l lw 2 lc "#D55E00" t "BDF2 (default)", \
        "timeSchemeRuns/newmark/postProcessing/0/solidPointDisplacement_pointDisp.dat" \
            u 1:5 w l lw 2 dt 2 lc "#0072B2" t "Newmark (trapezoidal)", \
        "timeSchemeRuns/bossak/postProcessing/0/solidPointDisplacement_pointDisp.dat" \
            u 1:5 w l lw 2 dt 3 lc "#009E73" t "Bossak-Newmark (alphaM = -0.1)"
} else {
    plot \
        "reference/abaqusC3D8.dat" u 1:2 w l lw 2 lc "black" t "Abaqus (C3D8)", \
        "postProcessing/0/solidPointDisplacement_pointDisp.dat" u 1:5 \
            w lp pt 6 ps 0.8 lw 2 lc "#D55E00" t sprintf("solids4foam (%s)", scheme)
}
