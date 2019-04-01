set xrange [0:0.005]

plot \
"energies.dat" u 1:2 w l t "external work",\
"" u 1:3 w l t "internal", \
"" u 1:4 w l t "kinetic", \
"" u 1:($3 + $4) w l t "internal + kinetic"