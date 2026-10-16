# Plot the modified large-deformation mesh study against Tukovic Fig. 28.
# The Allverify driver supplies data/output; defaults support manual use.
if (!exists("data")) data = "postProcessing/modified_mesh_sweep.csv"
if (!exists("output")) output = "postProcessing/modified_mesh_sweep.png"
baseDeltaX = 0.025

tukovicUx = 0.01463
tukovicUy = 0.005
tukovicUzSymmetryDifference = -0.000447

set datafile separator comma
set key left bottom opaque
set grid xtics ytics
set logscale x 2
set xrange [0.0025:0.03] reverse
set xtics ("0.003125" 0.003125, "0.00625" 0.00625, "0.00833" 0.00833333, "0.0125" 0.0125, "0.025" 0.025)
set format x "%.4g"
set xlabel "Solid mesh spacing {/Symbol D}x [m]"

set terminal pngcairo size 1800,900 enhanced font "Arial,16"
set output output
set multiplot layout 1,3 title "beamInCrossFlow modified case: mesh study at t = 8 s"
set ylabel "u_x(A) [m]"
plot data using (baseDeltaX*$6/0.05):9 with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    tukovicUx with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic Fig. 28"
set ylabel "u_y(A) [m]"
plot data using (baseDeltaX*$6/0.05):10 with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam", \
    tukovicUy with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic Fig. 28"
set format y "%.2e"
set ylabel "2 u_z(A) [m]"
plot data using (baseDeltaX*$6/0.05):15 with linespoints pt 7 ps 1.4 lw 2 lc rgb "#1f77b4" title "solids4foam (symmetry pair)", \
    tukovicUzSymmetryDifference with lines lw 2 dt 3 lc rgb "#6a3d9a" title "Tukovic Fig. 28"
unset multiplot
set output
