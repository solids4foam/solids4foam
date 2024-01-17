#!/usr/bin/gnuplot

set term pngcairo dashed font "Times,20" size 1920,1080

# Geometry input
R1 = 0.05
R2 = 0.07
R3 = 0.1

# Material input
E1 = 20e9
E2 = 200e9
Nu1 = 0.35
Nu2 = 0.3

# Internal pressure load
p = 1.0e5

set grid
set xrange [.050:0.100]
set xtics 0.01
set xlabel "r [m]"
set ylabel "sigma/p" font "Times Bold,20"
set key r t
set datafile separator " "
set style line 1 linecolor rgb 'red' linetype 6 linewidth 2 ps 2
set style line 2 linecolor rgb 'black' linetype 6 linewidth 2 ps 2
set sample 1000

pint = (2*R1**2*p/(E1*(R2**2-R1**2)))\
 / (((1.0/E2)*(((R3**2+R2**2)/(R3**2-R2**2))+Nu2))\
  + ((1.0/E1)*(((R2**2+R1**2)/(R2**2-R1**2))-Nu1)))

sigmatheta1(x) = (x < R2) ? (R1**2*p-R2**2*pint-(pint-p)*(R1*R2/x)**2)/(p*(R2**2-R1**2)) : 1/0
sigmatheta2(x) = (x > R2) ? (R2**2*pint+pint*(R2*R3/x)**2)/(p*(R3**2-R2**2)) : 1/0

sigmar1(x) = (x < R2) ? (R1**2*p-R2**2*pint+(pint-p)*(R1*R2/x)**2) / (p*(R2**2-R1**2)) : 1/0
sigmar2(x) = (x > R2) ? (R2**2*pint-pint*(R2*R3/x)**2) / (p*(R3**2-R2**2)) : 1/0

set output 'sigmaR.png'
set yrange [0:-1]
set ytics 0.1
set ylabel "{/Symbol s}_{r} / p_i" font "Times Bold,20"
plot sigmar1(x) w l ls 1 title"Analytical",\
     sigmar2(x) w l ls 1 notitle,\
    "postProcessing/sets/10/line_sigma:Transformed.xy" using 1:($2/p) w lp ls 2 notitle

set output 'sigmaTheta.png'
set yrange [-0.25:2.25]
set ytics 0.25
set ylabel "{/Symbol s}_{{/Symbol q}} / p_i" font "Times Bold,20"
plot sigmatheta1(x) w l ls 1 title"Analytical",\
     sigmatheta2(x) w l ls 1 notitle,\
     "postProcessing/sets/10/line_sigma\:Transformed.xy" using 1:($5/p) w lp ls 2 notitle
 

