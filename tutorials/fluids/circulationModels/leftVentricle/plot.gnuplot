# Plot the coupled left ventricle over one cardiac cycle (the last complete
# one by default, or the whole run if it is shorter than a cycle; another
# with: gnuplot -e "cycle=2" plot.gnuplot), and the volume exchanged between
# the ventricle and the circulation over the whole run.
#
# Writes:
#   interface.png            what the ventricle and the circulation exchange:
#                            the flow through the mitral (INLET) and aortic
#                            (OUTLET) valve planes, and the pressure the
#                            circulation imposes there with the left-atrial
#                            and aortic-sinus pressures each side of the valves
#   volumeConservation.png   the change in the ventricle volume against the
#                            change in the volume stored in the circulation,
#                            and their sum (volumeConservation.dat, written by
#                            Allrun)
#
# The results are in SI units (Pa, m3, m3/s) and are plotted in mmHg, mL and
# mL/s

dataFile = "postProcessing/circulationCoupling/0/coupling.dat"
volFile = "volumeConservation.dat"
haveVol = (system("test -s ".volFile." && echo 1 || echo 0") eq "1")

# Cardiac period (Tcycle in constant/circulationProperties)
Tcycle = 1.0

mmHg = 1.0/133.322387415
mL = 1e6

# Column of a quantity, from its name in the header line of the results file
col(name) = int(system(sprintf( \
    "awk -v n=%s '$1 == \"#\" && $2 == \"time\" " \
  . "{for (i = 3; i <= NF; i++) if ($i == n) {print i - 1; exit}}' %s", \
    name, dataFile)))

stats dataFile u 1 nooutput
if (!exists("cycle")) cycle = floor(STATS_max/Tcycle + 1e-6)
if (cycle >= 1) {
    t0 = (cycle - 1)*Tcycle
    t1 = cycle*Tcycle
    cycleTitle = sprintf("cycle %d", cycle)
} else {
    t0 = 0
    t1 = STATS_max
    cycleTitle = sprintf("t = 0 to %g s", t1)
}

# A quantity restricted to the cycle, scaled to the plotted units
inCycle(t) = (t >= t0 - 1e-9 && t <= t1 + 1e-9)
y(c, s) = inCycle($1) ? s*column(c) : NaN

set terminal pngcairo enhanced color font "Helvetica,11" size 800,700
set grid
set key top right font ",9"

# The ventricle-circulation interface
set output "interface.png"
set xrange [t0:t1]
set multiplot layout 2,1 title "Ventricle-circulation interface, ".cycleTitle
# The legends sit mid-cycle, where both valves are shut
set key at graph 0.43, graph 0.97 center top
set ylabel "Flow through the valve plane (mL/s)"
cMvQ = col("LV.inlet.Q")
cAvQ = col("LV.outlet.Q")
plot dataFile u 1:(y(cMvQ, mL)) w l lw 2 t "mitral, Q_{mv} (into the ventricle)", \
     dataFile u 1:(y(cAvQ, mL)) w l lw 2 t "aortic, Q_{av} (out of the ventricle)"
set xlabel "Time (s)"
set ylabel "Pressure (mmHg)"
cMvP = col("LV.inlet.Ppatch")
cAvP = col("LV.outlet.Ppatch")
cLAP = col("LA.P")
cSASP = col("SAS.C.P")
plot dataFile u 1:(y(cMvP, mmHg)) w l lw 2 lc 1 t "mitral plane, P_{mv}", \
     dataFile u 1:(y(cLAP, mmHg)) w l lw 1 dt 2 lc 1 t "left atrium", \
     dataFile u 1:(y(cAvP, mmHg)) w l lw 2 lc 2 t "aortic plane, P_{av}", \
     dataFile u 1:(y(cSASP, mmHg)) w l lw 1 dt 2 lc 2 t "aortic sinus"
unset multiplot

# Volume conservation over the whole run
if (haveVol) {
    set output "volumeConservation.png"
    set autoscale x
    set multiplot layout 2,1 title "Volume exchanged between the ventricle and the circulation"
    set key top right
    unset xlabel
    set ylabel "Volume change (mL)"
    plot volFile u 1:2 w l lw 4 lc rgb "black" t "ventricle, {/Symbol D}V_{LV}", \
         volFile u 1:3 w l lw 2 dt 2 lc rgb "red" t "circulation, -{/Symbol D}V_{0D}"
    set xlabel "Time (s)"
    set ylabel "{/Symbol D}V_{LV} + {/Symbol D}V_{0D} (mL)"
    set format y "%.1e"
    plot volFile u 1:4 w l lw 2 lc rgb "blue" notitle
    unset multiplot
}
