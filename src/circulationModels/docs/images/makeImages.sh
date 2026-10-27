#!/bin/bash
# Regenerate the figures of docs/TUTORIAL.md from the tutorial cases (needs gnuplot; the circuit diagrams need Python 3).
# Run the tutorial cases first (tutorials/tutorial/0*/ : circulationFoam; 04_cfdPipe: ./Allrun).
cd "$(dirname "$0")/../.." || exit 1
. /Volumes/OpenFOAM/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1 || true
T=tutorials/tutorial; I=docs/images
dat=postProcessing/circulation/0

export PATH=$PATH:$PWD/tools:/opt/homebrew/bin
for s in 01_windkessel 02_oneVentricleLoop; do
    (cd $T/$s && circulationFoam -diagram > /dev/null 2>&1)
    cp $T/$s/circuit.svg $I/step${s:1:1}_circuit.svg
done
(cd $T/04_cfdPipe/case && circulationFoam -diagram > /dev/null 2>&1) && cp $T/04_cfdPipe/case/circuit.svg $I/step4_circuit.svg

tools/circulationPlot -o $I/step1_pressure.png -t 4.2:5 $T/01_windkessel/$dat/circulation.dat Qin.out.P C.P
tools/circulationPlot -o $I/step2_pressures.png -t 19:20 $T/02_oneVentricleLoop/$dat/circulation.dat LV.P art.C.P ven.C.P
tools/circulationPlot -o $I/step2_flows.png -t 19:20 $T/02_oneVentricleLoop/$dat/circulation.dat AV.Q MV.Q

gnuplot -e "set terminal pngcairo size 700,550 font 'Helvetica,11'; set output '$I/step3_pvloops.png'; set termoption noenhanced; set grid; \
 set xlabel 'LV volume (mL)'; set ylabel 'LV pressure (mmHg)'; \
 plot '$T/02_oneVentricleLoop/$dat/LV.pvLoop.dat' using 1:2 with lines lw 2 title 'step 2: Emax 2.5 (one-ventricle loop)', \
      '$T/03b_presetOverride/$dat/LV.pvLoop.dat' using 1:2 with lines lw 2 title 'step 3b: healthyHuman, Emax_lv 1.5'"
