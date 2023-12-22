#!/usr/bin/gnuplot

set term pngcairo dashed font "Times, 28" size 1520,1080

set style line 1 linecolor rgb 'red' linetype 6 linewidth 2 ps 2
set style line 2 linecolor rgb 'black' linetype 6 linewidth 2 ps 2
set grid
set xrange [7:18.625]
set xlabel "r [m]"
set datafile separator " "
set sample 100

# Analytical solution
Ri = 7        # Inner cylinder radius
Ro = 18.625   # Outer cylinder radius
E = 1e10      # Young modulus
p = 10e6      # Internal pressure
nu = 0.495     # Poissons ratio

sigmaR(x) = (p*Ri**2 / (Ro**2-Ri**2))*(1-Ro**2/x**2)/p
sigmatheta(x) = (p*Ri**2 / (Ro**2-Ri**2))*(1+Ro**2/x**2)/p
dispR(x) = ((p/E) * (Ri**2/(Ro**2-Ri**2)))*((1-nu)*x+(1+nu)*(Ro**2/x))

# Path to postProcessing files
pathSigma="postProcessing/sets/10/line_sigma\:Transformed.xy"
pathD="postProcessing/sets/10/line_D.xy"

# Radial stress
set output 'sigmaR.png'
set ylabel "{/Symbol s}_{r} / p [-]"
set key r b
plot sigmaR(x) w l ls 1 title"Analytical",\
     pathSigma using 1:($2/p) w lp ls 2 title"solids4foam"

# Hoop stress
set output 'sigmaTheta.png'
set key r t
set ylabel "{/Symbol s}_{{/Symbol q}} / p [-]"
plot sigmatheta(x) w l ls 1 title"Analytical",\
     pathSigma using 1:($5/p) w lp ls 2 title"solids4foam"

# Radial displacement
set output 'dispR.png'
set key r t
set ylabel "u_r  [m]"
plot dispR(x) w l ls 1 title"Analytical",\
     pathD using 1:(sqrt($2**2+$3**2)) w lp ls 2 title"solids4foam"
