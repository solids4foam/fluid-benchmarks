# Plots the drag and the lift coefficients based on the data calculate by the
# `forceCoeffs` function object.

set terminal pdfcairo enhanced color solid

set datafile separator whitespace
set datafile commentschars "#"

# TODO(abzrg): Accept data file name from command-line
datafile = "./postProcessing/forceCoeffs/0/coefficient.dat" 
set output "forceCoeffsFnObject.pdf"

set title "Lift and Drag Coefficients vs. Time (via 'forceCoeffs' function object)"
set xlabel "Time, $t$ [s]"
set grid
set key top left

set ylabel "$C_d$ [1]"
set xrange [0:24]
set yrange [-4.5:4.5]
plot datafile using 1:2 title "Drag Coefficient" with lines

set ylabel "$C_l$ [1]"
set xrange [0:24]
set yrange [-0.06:0.06]
plot datafile using 1:3 title "Lift Coefficient" with lines
