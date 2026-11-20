# Plot one cardiac cycle of the coupled pipe closed loop (the last complete
# one by default, or the whole run if it is shorter than a cycle; another
# with: gnuplot -e "cycle=3" plot.gnuplot). The 0D reference of ./Allrun 0D
# is added to the ventricle plots when its results exist.
#
# Writes:
#   interface.png   what the pipe and the circuit exchange: the flow through
#                   the outlet and inlet patches and the pressure the circuit
#                   imposes there
#   heart.png       the ventricle: pressure-volume loop, chamber pressures,
#                   and the pressures each side of the aortic and mitral
#                   valves with the valve flows
#   circuit.png     the pressure in every compliance of the loop
#
# The results are in SI units (Pa, m3, m3/s) and are plotted in mmHg, mL and
# mL/s

dataFile = "postProcessing/circulationCoupling/0/coupling.dat"
refFile = "reference0D/postProcessing/circulation/0/circulation.dat"
haveRef = (system("test -f ".refFile." && echo 1 || echo 0") eq "1")

# Cardiac period (Tcycle in constant/circulationProperties)
Tcycle = 0.8

mmHg = 1.0/133.322387415
mL = 1e6

# Column of a quantity, from its name in the header line of a results file
colIn(name, file) = int(system(sprintf( \
    "awk -v n=%s '$1 == \"#\" && ($2 == \"time\" || $2 == \"t\") " \
  . "{for (i = 3; i <= NF; i++) if ($i == n) {print i - 1; exit}}' %s", \
    name, file)))
col(name) = colIn(name, dataFile)

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
set xrange [t0:t1]

# The pipe-circuit interface
set output "interface.png"
set multiplot layout 2,1 title "Pipe-circuit interface, ".cycleTitle
set ylabel "Flow through the patch (mL/s)"
cOutQ = col("outlet.Q")
cInQ = col("inlet.Q")
plot dataFile u 1:(y(cOutQ, mL)) w l lw 2 t "outlet (flow leaves the pipe)", \
     dataFile u 1:(y(cInQ, mL)) w l lw 2 dt 2 t "inlet (flow enters the pipe)"
set xlabel "Time (s)"
set ylabel "Pressure imposed by the circuit (mmHg)"
cOutP = col("outlet.Ppatch")
cInP = col("inlet.Ppatch")
plot dataFile u 1:(y(cOutP, mmHg)) w l lw 2 t "outlet (node C1)", \
     dataFile u 1:(y(cInP, mmHg)) w l lw 2 t "inlet (node C4)"
unset multiplot

# The heart
set terminal pngcairo enhanced color font "Helvetica,11" size 1100,800
set output "heart.png"
set multiplot layout 2,2 title "Heart, ".cycleTitle
cVP = col("ventricle.P")
cVV = col("ventricle.V")
cAP = col("atrium.P")
cC4P = col("C4.P")
if (haveRef) {
    rVP = colIn("ventricle.P", refFile)
    rVV = colIn("ventricle.V", refFile)
}

set title "Ventricle pressure-volume loop"
set autoscale x
set xlabel "Ventricle volume (mL)"
set ylabel "Ventricle pressure (mmHg)"
if (haveRef) {
    plot dataFile u (y(cVV, mL)):(y(cVP, mmHg)) w l lw 2 lc rgb "black" \
             t "coupled CFD", \
         refFile u (y(rVV, mL)):(y(rVP, mmHg)) w l lw 2 dt 2 \
             lc rgb "red" t "0D reference"
} else {
    plot dataFile u (y(cVV, mL)):(y(cVP, mmHg)) w l lw 2 lc rgb "black" \
        notitle
}

set title "Chamber pressures"
set xrange [t0:t1]
set xlabel "Time (s)"
set ylabel "Pressure (mmHg)"
if (haveRef) {
    plot dataFile u 1:(y(cAP, mmHg)) w l lw 2 t "atrium", \
         dataFile u 1:(y(cVP, mmHg)) w l lw 2 t "ventricle", \
         dataFile u 1:(y(cC4P, mmHg)) w l lw 2 t "aortic side (C4)", \
         refFile u 1:(y(rVP, mmHg)) w l lw 1 dt 2 lc rgb "red" \
             t "ventricle, 0D reference"
} else {
    plot dataFile u 1:(y(cAP, mmHg)) w l lw 2 t "atrium", \
         dataFile u 1:(y(cVP, mmHg)) w l lw 2 t "ventricle", \
         dataFile u 1:(y(cC4P, mmHg)) w l lw 2 t "aortic side (C4)"
}

set y2label "Valve flow (mL/s)"
set ytics nomirror
set y2tics
cAoQ = col("aortic.Q")
set title "Aortic valve: ventricle | C4"
plot dataFile u 1:(y(cVP, mmHg)) w l lw 2 t "ventricle (upstream)", \
     dataFile u 1:(y(cC4P, mmHg)) w l lw 2 t "C4 (downstream)", \
     dataFile u 1:(y(cAoQ, mL)) axes x1y2 w l dt 2 lc rgb "gray40" \
         t "aortic.Q"
cMiQ = col("mitral.Q")
set title "Mitral valve: atrium | ventricle"
plot dataFile u 1:(y(cAP, mmHg)) w l lw 2 t "atrium (upstream)", \
     dataFile u 1:(y(cVP, mmHg)) w l lw 2 t "ventricle (downstream)", \
     dataFile u 1:(y(cMiQ, mL)) axes x1y2 w l dt 2 lc rgb "gray40" \
         t "mitral.Q"
unset multiplot
unset title
unset y2label
unset y2tics
set ytics mirror

# The pressure in every compliance
set terminal pngcairo enhanced color font "Helvetica,11" size 800,500
set output "circuit.png"
set title "Compliance pressures, ".cycleTitle
set xlabel "Time (s)"
set ylabel "Pressure (mmHg)"
cC1P = col("C1.P")
cC2P = col("C2.P")
cC3P = col("C3.P")
plot dataFile u 1:(y(cC1P, mmHg)) w l lw 2 t "C1 (pipe outlet)", \
     dataFile u 1:(y(cC2P, mmHg)) w l lw 2 t "C2", \
     dataFile u 1:(y(cC3P, mmHg)) w l lw 2 t "C3", \
     dataFile u 1:(y(cC4P, mmHg)) w l lw 2 t "C4 (pipe inlet)"

# Ventricle volumes and peak pressure of the cycle
stats dataFile u (y(cVV, mL)) nooutput
edv = STATS_max
esv = STATS_min
stats dataFile u (y(cVP, mmHg)) nooutput
set print "-"
print sprintf( \
    "%s: ventricle EDV %.1f mL, ESV %.1f mL, SV %.1f mL, " \
  . "peak ventricle pressure %.0f mmHg", \
    cycleTitle, edv, esv, edv - esv, STATS_max)
