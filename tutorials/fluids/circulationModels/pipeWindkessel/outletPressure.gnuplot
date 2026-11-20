# Outlet pressure of the coupled run against the exact solution
#   Pc(t) = A R/(1 + (wRC)^2) [sin wt - wRC cos wt + wRC exp(-t/RC)]
#   Po(t) = Pc(t) + Rch A sin wt
# with A = 1, w = 2 pi, R = 10, Rch = 1, C = 0.01 and Pc(0) = 0
set terminal pngcairo enhanced color font "Helvetica,11" size 800,500
set output "outletPressure.png"

dataFile = "postProcessing/circulationCoupling/0/coupling.dat"

A = 1.0
w = 2*pi
R = 10.0
Rch = 1.0
C = 0.01
k = w*R*C
Pc(t) = A*R/(1 + k**2)*(sin(w*t) - k*cos(w*t) + k*exp(-t/(R*C)))
Po(t) = Pc(t) + Rch*A*sin(w*t)

set xlabel "Time (s)"
set ylabel "Outlet pressure P_o (Pa)"
set grid
set key top right
set xrange [0:1]
set samples 500

# Column 4 is outlet.Ppatch
plot Po(x) w l lw 2 lc rgb "black" t "Exact", \
     dataFile u 1:4 w p pt 7 ps 1 lc rgb "red" t "solids4Foam"
