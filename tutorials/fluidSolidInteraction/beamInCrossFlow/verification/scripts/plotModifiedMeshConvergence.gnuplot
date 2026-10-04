# Plot the modified large-deformation mesh study against Tukovic Fig. 28.
# The Allverify driver supplies data/output; defaults support manual use.
if (!exists("data")) data = "postProcessing/modified_mesh_sweep.csv"
if (!exists("output")) output = "postProcessing/modified_mesh_sweep.png"
tukovicUx = 0.01463
tukovicUy = 0.005
tukovicUzSymmetryDifference = -0.000447

set datafile separator comma
set key left bottom opaque
set grid xtics ytics
set logscale x 2
set format x "%.4g"
set xlabel "Representative near-body cell size [m]"

set terminal pngcairo size 1800,900 enhanced font "Arial,16"
set output output
set multiplot layout 1,3 title "beamInCrossFlow modified case: mesh study at t = 8 s"
set ylabel "u_x(A) [m]"
plot data using (column("near_body_cell_size")):(column("ux")) with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    tukovicUx with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic Fig. 28"
set ylabel "u_y(A) [m]"
plot data using (column("near_body_cell_size")):(column("uy")) with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    tukovicUy with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic Fig. 28"
set format y "%.2e"
set ylabel "2 u_z(A) [m]"
plot data using (column("near_body_cell_size")):(column("uz_symmetry_difference")) with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam (symmetry pair)", \
    tukovicUzSymmetryDifference with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic Fig. 28"
unset multiplot
set output
