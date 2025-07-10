set term pdfcairo dashed enhanced
set datafile separator " "

if (ARGC < 1) {
    print "Error: No input configuration name provided."
    print "usage: ", ARG0, " <configName>"
    exit
} else {
    # Summary of all configurations are combined into one data file
    datafile = ARG1
}

set output "linearFitVelocityErrors.pdf"
slopfile = "orderOfAccuracySlopesVelocity.dat"

set grid
set xtics add (25, 50, 100, 200)
set ytics
set xrange [10:200]
set yrange [1e-5:0.1]
set logscale x
set logscale y
set format y "10^{%L}"
set xlabel "Average cell spacing (in mm)"
set ylabel "Error (in m/s)"
set key right bottom;
set key spacing 2 # Adds space between legend entries

set label "1^{st} order" at graph 0.43,0.76 center rotate by 10
set label "2^{nd} order" at graph 0.43,0.33 center rotate by 22

# Average mesh spacing of mesh1
dx = 0.2

# A linear function for fitting in log-log space
g(x) = m * x + c

# Silence verbose fit output
set fit quiet

# Fit g(x) to the L_1 error data
set fit logfile "L_1.fit.log"
fit [0:1000] g(x) datafile using (log(1e3*dx/(2**($1 - 1)))):(log($4)) via m, c
plot datafile using (1e3*dx/(2**($1 - 1))):4 with points pt 6 lc "red" title "L_1", \
     x**m * exp(c) with lines lw 2 lc rgb "dark-red" title sprintf("Y = {%.2f}X + %.2g", m, c), \
     slopfile u 1:2 w l lw 2 lc "gray" notitle, \
     slopfile u 1:3 w l lw 2 lc "gray" notitle

# Fit g(x) to the L_2 error data
set fit logfile "L_2.fit.log"
fit [0:1000] g(x) datafile using (log(1e3*dx/(2**($1 - 1)))):(log($5)) via m, c
plot datafile using (1e3*dx/(2**($1 - 1))):5 with points pt 6 lc "green" title "L_2", \
     x**m * exp(c) with lines lw 2 lc rgb "dark-green" title sprintf("Y = {%.2f}X + %.2g", m, c), \
     slopfile u 1:2 w l lw 2 lc "gray" notitle, \
     slopfile u 1:3 w l lw 2 lc "gray" notitle

# Fit g(x) to the L_∞ error data
set fit logfile "L_inf.fit.log"
fit [0:1000] g(x) datafile using (log(1e3*dx/(2**($1 - 1)))):(log($6)) via m, c
plot datafile using (1e3*dx/(2**($1 - 1))):6 with points pt 6 lc "blue" title "L_∞", \
     x**m * exp(c) with lines lw 2 lc rgb "dark-blue" title sprintf("Y = {%.2f}X + %.2g", m, c), \
     slopfile u 1:2 w l lw 2 lc "gray" notitle, \
     slopfile u 1:3 w l lw 2 lc "gray" notitle
