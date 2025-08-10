# Plots the drag and the lift coefficients based on the data calculate by the
# `forceCoeffs` function object.

set terminal pdfcairo enhanced color solid

set datafile commentschars "#"
set datafile separator whitespace
dataFile = "./postProcessing/forceCoeffs/0/coefficient.dat"
validCdDataFile = "./Cd.dat"
validClDataFile = "./Cl.dat"

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

set title "Drag Coefficient (C_d) vs. Time (t)" . \
          " (using {/Monospace forceCoeffs} function object)"

set xlabel "t [s]"
set ylabel "C_d [1]"

set xtics 2
set ytics 1.5

set xrange [0:24]
set yrange [-4.5:4.5]

plot \
    dataFile \
        using 1:2 \
        title "Present" \
        with lines lc rgb lineColor, \
    validCdDataFile \
        using 1:2 \
        title "Wan\\&Turek" \
        with points pt pointType ps pointSize lc rgb pointColor

# --- Plot Lift Coefficient ---------------------------------------------------#

set title "Lift Coefficient (C_l) vs. Time (t)" . \
          " (using the {/Monospace forceCoeffs} function object)"

set xlabel "t [s]"
set ylabel "C_l [1]"

set xtics 2
set ytics 0.02

set xrange [0:24]
set yrange [-0.06:0.06]

plot \
    dataFile \
        using 1:3 \
        title "Present" \
        with lines lc rgb lineColor, \
    validClDataFile \
        using 1:2 \
        title "Wan\\&Turek" \
        with points pt pointType ps pointSize lc rgb pointColor
