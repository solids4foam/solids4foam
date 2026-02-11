reset
#set terminal pdfcairo enhanced dashed
set terminal pdfcairo enhanced dashed size 20cm,20cm
set output "displacement.pdf"

set multiplot

# ==================================================
# 1) Background image (full page, untouched)
# ==================================================
unset tics
unset border
unset xlabel
unset ylabel
unset key

img_ratio = 945.0 / 1196.0

bg_w = 0.9
bg_h = bg_w * img_ratio

bg_x0 = (1 - bg_w) / 2.0
bg_y0 = (1 - bg_h) / 2.0

set origin bg_x0, bg_y0
set size   bg_w,  bg_h

plot "images/Cerquaglia-fillingElasticContainer-apex-displacement.png" \
     binary filetype=png with rgbimage notitle

# ==================================================
# 2) Foreground plot:
# ==================================================

# Position and scale parameters for aligning the plot borders
fg_x0 = 0.169   # left
fg_y0 = 0.248   # bottom
fg_w  = 0.755   # width
fg_h  = 0.59   # height
# -------------------------------------------

set origin fg_x0, fg_y0
set size   fg_w,  fg_h

set xrange [0:10]
set yrange [-1.2:0.8]

unset grid
unset key

set border

set key reverse
set key Left
set key font ",18"
set key at screen 0.95,0.68
set key spacing 1.2

plot \
  "postProcessing/0/solidPointDisplacement_disp.dat" \
    u 1:3 w l lw 4 dt 2 lc "red" t "solids4foam (480+2000 cells)"

unset multiplot
set output