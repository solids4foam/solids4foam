# Plot beamInCrossFlow mesh-verification predictions against references.
# The Allverify driver supplies data/output; defaults support manual use.
if (!exists("data")) data = "postProcessing/original_mesh_sweep.csv"
if (!exists("output")) output = "postProcessing/original_mesh_sweep.png"
baseDeltaX = 0.025

richterUx = 5.95e-5
tukovicUx = 5.93e-5
tukovicUy = 2.40e-5
richterFx = 1.33
tukovicFx = 1.31
tukovicFy = 0.1055

set datafile separator comma
set key left bottom opaque
set grid xtics ytics
set logscale x 2
set xrange [0.0025:0.03] reverse
set xtics ("0.003125" 0.003125, "0.00625" 0.00625, "0.0125" 0.0125, "0.025" 0.025)
set format x "%.4g"
set xlabel "Solid mesh spacing {/Symbol D}x [m]"

set terminal pngcairo size 1800,1200 enhanced font "Arial,16"
set output output
set multiplot layout 2,2 title "beamInCrossFlow original case: mesh study at t = 8 s"
set format y "%.1e"
set ylabel "u_x(A) [m]"
plot data using (baseDeltaX/(2.0**$4)):9 with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    richterUx with lines lw 2 dt 2 lc rgb "#cc0000" title "Richter benchmark", \
    tukovicUx with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic OpenFOAM"
set ylabel "u_y(A) [m]"
plot data using (baseDeltaX/(2.0**$4)):10 with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    tukovicUy with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic OpenFOAM"
set ylabel "F_x [N]"
set format y "%.3f"
plot data using (baseDeltaX/(2.0**$4)):12 with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    richterFx with lines lw 2 dt 2 lc rgb "#cc0000" title "Richter benchmark", \
    tukovicFx with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic OpenFOAM"
set ylabel "F_y [N]"
plot data using (baseDeltaX/(2.0**$4)):13 with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    tukovicFy with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic OpenFOAM"
unset multiplot
set output
