# Plot the tip displacement history against the Abaqus reference solution
set term pngcairo dashed size 1024,768 font "Arial,18"
set output "tipDisplacement.png"

set grid
set xlabel "Time (s)"
set ylabel "Tip displacement magnitude (m)"
set key right top

# Show the reference over one full period, which matches the tutorial end time
set xrange [0:0.65]

plot \
    "reference/abaqusC3D8.dat" u 1:2 w l lw 2 lc "black" t "Abaqus (C3D8)", \
    "postProcessing/0/solidPointDisplacement_pointDisp.dat" u 1:5 \
        w lp pt 6 ps 0.8 lw 2 lc "red" t "solids4foam"
