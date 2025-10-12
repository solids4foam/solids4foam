set xlabel "Time (in s)"
set ylabel "Number of thermo-FSI iterations per time-step"
set grid
while (1){
    plot "residuals/thermalResiduals.dat" u 1:2 w l
    pause 1
}