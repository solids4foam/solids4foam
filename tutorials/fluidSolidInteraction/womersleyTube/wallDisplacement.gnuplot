set terminal pdfcairo enhanced color solid font "Helvetica,10" size 5,3.5
set output "wallDisplacement.pdf"
set xlabel "Time (s)"
set ylabel "Radial wall displacement (mm)"
set grid
set key outside top center horizontal

# Exact solution, Re[eta exp(i (omega t - k x))], at x = L/4, L/2 and 3L/4
# (see verification/reference/womersleyTube_verification_references.json)
omega = 2*pi*0.02
exact(t, a, phi) = 1e3*a*cos(omega*t + phi)

# The points are on the wedge plane at -0.5 degrees
c = cos(0.5*pi/180)
s = sin(0.5*pi/180)
radial(dy, dz) = 1e3*(dy*c - dz*s)

plot \
    "postProcessing/0/solidPointDisplacement_wallQuarter.dat" \
        u 1:(radial($3, $4)) w l lw 2 lc rgb "#1b9e77" t "x = L/4", \
    exact(x, 4.877217512946812e-4, -0.5053998829035038) \
        w l dt 2 lw 1.5 lc rgb "black" t "Exact", \
    "postProcessing/0/solidPointDisplacement_wallMid.dat" \
        u 1:(radial($3, $4)) w l lw 2 lc rgb "#d95f02" t "x = L/2", \
    exact(x, 4.512494768325512e-4, -1.0453187730902433) \
        w l dt 2 lw 1.5 lc rgb "black" notitle, \
    "postProcessing/0/solidPointDisplacement_wallThreeQuarter.dat" \
        u 1:(radial($3, $4)) w l lw 2 lc rgb "#7570b3" t "x = 3L/4", \
    exact(x, 4.175046321004123e-4, -1.585237663276983) \
        w l dt 2 lw 1.5 lc rgb "black" notitle
