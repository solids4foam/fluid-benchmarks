set term pdfcairo dashed enhanced
set datafile separator " "

set output "velocityErrors.pdf"

#set size ratio 1

set grid
set xrange [10:200]
set yrange [1e-5:0.1]
set xtics add (25, 50, 100, 200)
set ytics
set logscale x
set logscale y
set format y "10^{%L}"
#set ytics 0.002
set xlabel "Average cell spacing (in mm)"
set ylabel "Error (in m/s)"
set key right bottom;

set label "1^{st} order" at graph 0.5,0.78 center rotate by 10
set label "2^{nd} order" at graph 0.5,0.37 center rotate by 25

# Average mesh spacing of mesh1
dx=0.2

# 15,14 - upturned solid penta
# 5,4 - square

# Assume the mesh spacing is being halved for each succesive mesh
plot \
    "hex.lu.summary.txt" u (1e3*dx/(2**($0))):($4) w lp pt 5 lc "green" t "L_1", \
    "hex.lu.summary.txt" u (1e3*dx/(2**($0))):($5) w lp pt 15 lc "red" t "L_2", \
    "hex.lu.summary.txt" u (1e3*dx/(2**($0))):($6) w lp pt 9 lc "blue" t "L_∞", \
    "orderOfAccuracySlopesVelocity.dat" u 1:2 w l lw 2 lc "black" notitle, \
    "orderOfAccuracySlopesVelocity.dat" u 1:3 w l lw 2 lc "black" notitle
