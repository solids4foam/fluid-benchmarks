set term pdfcairo dashed enhanced
set datafile separator " "

set output "velocity_orderOfAccuracy.pdf"

set grid
set xrange [10:200]
set yrange [0:3]
set xtics
set xtics add (25, 50, 100, 200)
set ytics
set logscale x
#set logscale y
#set ytics 0.002
set xlabel "Average cell spacing (in mm)"
set ylabel "Order of accuracy"
set key bottom left;

# Average mesh spacing of mesh1
dx=0.2

# Assume the mesh spacing is being halved for each succesive mesh
plot \
    "hex.lu.orderOfAccuracy.txt" using (1e3*dx/(2**($0))):2 skip 1 w lp pt 5 lc "green" t "L_1", \
    "hex.lu.orderOfAccuracy.txt" using (1e3*dx/(2**($0))):3 skip 1 w lp pt 5 lc "red" t "L_2", \
    "hex.lu.orderOfAccuracy.txt" using (1e3*dx/(2**($0))):4 skip 1 w lp pt 4 lc "blue" t "L_∞"

