#!/usr/bin/gnuplot

set term pngcairo dashed font "Times,20" size 1920,1080

# Thickness of the mid-layer relative to W: Allrun sets it with
# gnuplot -e "hOverW=<value>"
if (!exists("hOverW")) hOverW = 0.01

# Geometry input
W = 0.03

# Material input
E1 = 206e9
E2 = 70.3e9
E3 = 4.93e9
Nu1 = 0.3
Nu2 = 0.345
Nu3 = 0.33

# Tensile load
sigma0 = 1.0e6

# Order p of the stress singularity, sigma ~ r^(-p), where the interface
# between materials a and b meets the free edge, in plane strain. It is the
# root lambda = 1 - p of the characteristic equation of Bogy, written with the
# Dundurs' parameters alpha and beta, which is found by bisection between
# lambda = 0.5 and the trivial root lambda = 1
mu(E, nu) = E/(2*(1 + nu))
m(nu) = 4*(1 - nu)
k(Ea, nua, Eb, nub) = mu(Eb, nub)/mu(Ea, nua)
alpha(Ea, nua, Eb, nub) = (m(nub) - k(Ea, nua, Eb, nub)*m(nua))\
 / (m(nub) + k(Ea, nua, Eb, nub)*m(nua))
beta(Ea, nua, Eb, nub) = ((m(nub) - 2) - k(Ea, nua, Eb, nub)*(m(nua) - 2))\
 / (m(nub) + k(Ea, nua, Eb, nub)*m(nua))
s(l) = sin(pi*l/2)**2
f(l, a, b) = l**2*(l**2 - 1)*a**2 + 2*l**2*(s(l) - l**2)*a*b\
 + (s(l) - l**2)**2*b**2 + s(l)*(1 - s(l))
root(lo, hi, a, b, n) = (n == 0) ? 0.5*(lo + hi)\
 : (f(lo, a, b)*f(0.5*(lo + hi), a, b) <= 0\
 ? root(lo, 0.5*(lo + hi), a, b, n - 1)\
 : root(0.5*(lo + hi), hi, a, b, n - 1))
p(Ea, nua, Eb, nub) = 1 - root(0.5, 0.9999, alpha(Ea, nua, Eb, nub),\
 beta(Ea, nua, Eb, nub), 40)

p12 = p(E1, Nu1, E2, Nu2)
p13 = p(E1, Nu1, E3, Nu3)
p23 = p(E2, Nu2, E3, Nu3)

set grid
set logscale xy
set xrange [1e-5:1]
set format x "10^{%L}"
set xlabel "r / W"
set ylabel "{/Symbol s}_{xx} / {/Symbol s}_0" font "Times Bold,20"
set key r t
set datafile separator whitespace
set style line 1 linecolor rgb 'red' linetype 1 linewidth 2
set style line 2 linecolor rgb 'black' linetype 6 linewidth 2 ps 2
set style line 3 linecolor rgb 'blue' linetype 1 linewidth 2
set style line 4 linecolor rgb 'black' linetype 4 linewidth 2 ps 2

# The singular term is drawn with the analytical order and an amplitude fitted
# to the numerical stress in rMin < r/W < rMax, which is close enough to the
# free edge for the singular term to dominate, but more than a few cells
# from it
rMin = 1e-4
rMax = 3e-4
set fit quiet
set fit logfile "/dev/null"

set output 'sigmaxx.png'

if (hOverW > 0) {
    file13 = "postProcessing/sample/1/interface13_sigma.xy"
    file23 = "postProcessing/sample/1/interface23_sigma.xy"

    lnK13 = 1
    lnK23 = 1
    fit [rMin:rMax] lnK13 - p13*log(x) file13\
        using (($1 + W)/W):(log($2/sigma0)) via lnK13
    fit [rMin:rMax] lnK23 - p23*log(x) file23\
        using (($1 + W)/W):(log($2/sigma0)) via lnK23

    set yrange [0.7:4]
    set ytics (0.7, 0.8, 0.9, 1, 1.5, 2, 3, 4)
    plot exp(lnK13)*x**(-p13) w l ls 1\
            title sprintf("Analytical order, p_{13} = %.4f", p13),\
         file13 using (($1 + W)/W):($2/sigma0) w p ls 2\
            title "Interface 13",\
         exp(lnK23)*x**(-p23) w l ls 3\
            title sprintf("Analytical order, p_{23} = %.4f", p23),\
         file23 using (($1 + W)/W):($2/sigma0) w p ls 4\
            title "Interface 23"
} else {
    file12 = "postProcessing/sample/1/interface12_sigma.xy"

    lnK12 = 1
    fit [rMin:rMax] lnK12 - p12*log(x) file12\
        using (($1 + W)/W):(log($2/sigma0)) via lnK12

    set yrange [0.9:3]
    set ytics (0.9, 1, 1.5, 2, 3)
    plot exp(lnK12)*x**(-p12) w l ls 1\
            title sprintf("Analytical order, p_{12} = %.4f", p12),\
         file12 using (($1 + W)/W):($2/sigma0) w p ls 2\
            title "Interface 12"
}
