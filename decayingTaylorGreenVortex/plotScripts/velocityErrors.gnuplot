set term pdfcairo dashed enhanced
set datafile separator " "

if (ARGC < 1) {
    print "Error: No input configuration name provided."
    print "usage: ", ARG0, " <configName>"
    exit
} else {
    configName = ARG1
}

set output configName.".velocityErrors.pdf"

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
dx = 0.2

# 15,14 - upturned solid penta
# 5,4 - square

# Assume the mesh spacing is being halved for each succesive mesh
datafile = configName.".summary.txt"
slopfile = "orderOfAccuracySlopesVelocity.dat"
plot \
    datafile u (1e3*dx/(2**($0))):($4) w lp pt 5 lc "green" t "L_1", \
    datafile u (1e3*dx/(2**($0))):($5) w lp pt 15 lc "red" t "L_2", \
    datafile u (1e3*dx/(2**($0))):($6) w lp pt 9 lc "blue" t "L_∞", \
    slopfile u 1:2 w l lw 2 lc "black" notitle, \
    slopfile u 1:3 w l lw 2 lc "black" notitle
