# Flow rate through the valve (the outlet), in ml/s, with the phases of the
# valve motion: open, closing, closed and opening
set terminal pdfcairo enhanced color font "Helvetica,12" size 6,4
set output "flowRate.pdf"

flowFile = "postProcessing/flowRate/0/surfaceFieldValue.dat"

set xlabel "Time (s)"
set ylabel "Flow rate (ml/s)"
set grid
set key top right

# Closing and opening windows of the sliceAxis valve motion (period 0.8 s)
set object 1 rect from 0.24, graph 0 to 0.32, graph 1 fc rgb "#eeeeee" behind
set object 2 rect from 0.64, graph 0 to 0.76, graph 1 fc rgb "#eeeeee" behind
set object 3 rect from 1.04, graph 0 to 1.12, graph 1 fc rgb "#eeeeee" behind
set object 4 rect from 1.44, graph 0 to 1.56, graph 1 fc rgb "#eeeeee" behind

plot flowFile u 1:(1e6*$2) w l lw 2 lc rgb "black" t "Outlet"
