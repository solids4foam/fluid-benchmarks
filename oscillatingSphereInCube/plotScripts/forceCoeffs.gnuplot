# Plots the drag and the lift coefficients based on the data calculate by the
# `forceCoeffs` function object.

set terminal pdfcairo enhanced color solid

set datafile separator whitespace
set datafile commentschars "#"
dataFile = "./postProcessing/forceCoeffs/0/coefficient.dat"
CdDataFile = "./Cd.dat"

set output "forceCoeffsFnObject.pdf"

set grid
set tics nomirror
set border 3
set key top left
set key box width 2 height 1

lineColor = "blue"
pointColor = "red"
pointType = 6 # Hollow circle
pointSize = .2

# --- Plot Drag Coefficient ---------------------------------------------------#

set title "Drag Coefficient (C_d) vs. Time (t)"

set xlabel "t [s]"
set ylabel "C_d [1]"

set xtics 1
set ytics 10

set xrange [0:10]
set yrange [-30:30]

plot \
    dataFile \
        using 1:2 \
        title "Present" \
        with lines lc rgb lineColor, \
    CdDataFile \
        using 1:2 \
        title "Erzincanli\\&Sahin" \
        with points pt pointType ps pointSize lc rgb pointColor
