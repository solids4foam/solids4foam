# Plot one cardiac cycle of the healthy human circulation (the last one by
# default; another with: gnuplot -e "cycle=12" plot.gnuplot)
#
# Writes:
#   pvLoops.png         pressure-volume loops of the four chambers
#   valvePressures.png  the pressure each side of the four valves, with the
#                       valve flow
#   wiggers.png         the left heart: pressures, left ventricular volume and
#                       the aortic and mitral flows
#
# The results are in mmHg, mL and mL/s (units mmHg_mL in
# constant/circulationProperties)

dataFile = "postProcessing/circulation/0/circulation.dat"

# Cardiac period (Tcycle in constant/circulationProperties)
Tcycle = 1.0

# Column of a quantity, from its name in the "# t ..." header line
col(name) = int(system(sprintf( \
    "awk -v n=%s '$1 == \"#\" && $2 == \"t\" " \
  . "{for (i = 3; i <= NF; i++) if ($i == n) {print i - 1; exit}}' %s", \
    name, dataFile)))

stats dataFile u 1 nooutput
if (!exists("cycle")) cycle = floor(STATS_max/Tcycle + 1e-6)
t0 = (cycle - 1)*Tcycle
t1 = cycle*Tcycle

# Time in the cycle, and a quantity restricted to the cycle
inCycle(t) = (t >= t0 - 1e-9 && t <= t1 + 1e-9)
tc = '(inCycle($1) ? $1 - t0 : NaN)'
y(c) = inCycle($1) ? column(c) : NaN

set terminal pngcairo enhanced color font "Helvetica,11" size 1100,800
set grid
set key top right font ",9"

# Pressure-volume loops
set output "pvLoops.png"
set multiplot layout 2,2 title sprintf("Pressure-volume loops, cycle %d", cycle)
do for [ch in "LV RV LA RA"] {
    cV = col(ch.".V")
    cP = col(ch.".P")
    set title ch
    set xlabel "Volume (mL)"
    set ylabel "Pressure (mmHg)"
    plot dataFile u (y(cV)):(y(cP)) w l lw 2 lc rgb "black" notitle
}
unset multiplot

# The pressure each side of the valves, with the valve flow
set output "valvePressures.png"
set multiplot layout 2,2 \
    title sprintf("Pressure each side of the valves, cycle %d", cycle)
set xlabel "Time in cycle (s)"
set ylabel "Pressure (mmHg)"
set y2label "Valve flow (mL/s)"
set ytics nomirror
set y2tics
upList = "LV.P LA.P RV.P RA.P"
downList = "SAS.C.P LV.P PAS.C.P RV.P"
valveList = "AV MV PV TV"
titleList = "'Aortic valve: LV | aorta (SAS)' \
'Mitral valve: LA | LV' \
'Pulmonary valve: RV | pulmonary artery (PAS)' \
'Tricuspid valve: RA | RV'"
do for [k = 1:4] {
    up = word(upList, k)
    down = word(downList, k)
    cUp = col(up)
    cDown = col(down)
    cQ = col(word(valveList, k).".Q")
    set title word(titleList, k)
    plot dataFile u @tc:(y(cUp)) w l lw 2 t up." (upstream)", \
         dataFile u @tc:(y(cDown)) w l lw 2 t down." (downstream)", \
         dataFile u @tc:(y(cQ)) axes x1y2 w l dt 2 lc rgb "gray40" \
             t word(valveList, k).".Q"
}
unset multiplot
unset y2label
unset y2tics
set ytics mirror

# Wiggers diagram of the left heart
set terminal pngcairo enhanced color font "Helvetica,11" size 800,900
set output "wiggers.png"
set multiplot layout 3,1 title sprintf("Left heart, cycle %d", cycle)
unset title
unset xlabel
set ylabel "Pressure (mmHg)"
cLVP = col("LV.P")
cAoP = col("SAS.C.P")
cLAP = col("LA.P")
plot dataFile u @tc:(y(cLVP)) w l lw 2 t "LV", \
     dataFile u @tc:(y(cAoP)) w l lw 2 t "aorta (SAS)", \
     dataFile u @tc:(y(cLAP)) w l lw 2 t "LA"
set ylabel "LV volume (mL)"
cLVV = col("LV.V")
plot dataFile u @tc:(y(cLVV)) w l lw 2 lc rgb "forest-green" notitle
set xlabel "Time in cycle (s)"
set ylabel "Flow (mL/s)"
cAVQ = col("AV.Q")
cMVQ = col("MV.Q")
plot dataFile u @tc:(y(cAVQ)) w l lw 2 lc rgb "red" t "aortic valve", \
     dataFile u @tc:(y(cMVQ)) w l lw 2 lc rgb "purple" t "mitral valve"
unset multiplot

# Left ventricular volumes of the cycle
stats dataFile u (y(cLVV)) nooutput
set print "-"
print sprintf( \
    "cycle %d: LV EDV %.1f mL, ESV %.1f mL, SV %.1f mL, EF %.1f %%", \
    cycle, STATS_max, STATS_min, STATS_max - STATS_min, \
    100*(STATS_max - STATS_min)/STATS_max)
