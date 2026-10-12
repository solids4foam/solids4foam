# Plot beamInCrossFlow mesh-verification predictions against references.
# The Allverify driver supplies data/output; defaults support manual use.
if (!exists("data")) data = "postProcessing/richter_graded_mesh_sweep.csv"
if (!exists("output")) output = "postProcessing/richter_graded_mesh_sweep.png"
richterUx = 5.924e-5
richterFx = 1.327

set datafile separator comma
set key left bottom opaque
set grid xtics ytics
set logscale x 2
set format x "%.4g"
set xlabel "Representative near-body cell size [m]"

set terminal pngcairo size 1800,1200 enhanced font "Arial,16"
set output output
set multiplot layout 2,2 title "beamInCrossFlow Richter benchmark: graded-mesh study"
set format y "%.1e"
set ylabel "u_x(A) [m]"
plot data using (column("near_body_cell_size")):(column("ux")) with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    richterUx with lines lw 2 dt 2 lc rgb "#cc0000" title "Richter Table 6"
set ylabel "u_y(A) [m]"
plot data using (column("near_body_cell_size")):(column("uy")) with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam"
set ylabel "F_x [N]"
set format y "%.3f"
plot data using (column("near_body_cell_size")):(column("fx")) with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    richterFx with lines lw 2 dt 2 lc rgb "#cc0000" title "Richter Table 6"
set ylabel "F_y [N]"
plot data using (column("near_body_cell_size")):(column("fy")) with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam"
unset multiplot
set output
