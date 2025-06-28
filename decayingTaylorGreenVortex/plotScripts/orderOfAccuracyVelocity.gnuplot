set term pdfcairo dashed enhanced
set datafile separator " "

if (ARGC < 1) {
    print "Error: No input configuration name provided."
    print "usage: ", ARG0, " <configName>"
    exit
} else {
    configName = ARG1
}

set output configName.".velocity_orderOfAccuracy.pdf"

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
    configName.".orderOfAccuracy.txt" using (1e3*dx/(2**($0))):2 skip 1 w lp pt 5 lc "green" t "L_1", \
    configName.".orderOfAccuracy.txt" using (1e3*dx/(2**($0))):3 skip 1 w lp pt 5 lc "red" t "L_2", \
    configName.".orderOfAccuracy.txt" using (1e3*dx/(2**($0))):4 skip 1 w lp pt 4 lc "blue" t "L_∞"

