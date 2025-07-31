# Plots the drag and the lift coefficients based on the data calculate by the
# `forces` function object.

set terminal pdfcairo enhanced color solid

set datafile separator whitespace
set datafile commentschars "#"

# TODO(abzrg): Accept data file name from command-line
datafile =  "./postProcessing/forces/0/force.dat"
set output 'forceCoeffsForces.pdf'

set title "Lift and Drag Coefficients vs. Time (via 'forces' function object)"
set xlabel "Time, $t$ [s]"
set grid
set key top left

# ------------------------------- Parameters -------------------------------- #

# Frequency of oscillation of displacement of cylinder
f = 0.25
# Magnitude of oscillation of displacement of cylinder
A = 0.25
# Magnitude of reference velocity
magUInf = 2*pi*f*A

# Span of mesh in the z (depth) direction
# See blockMeshDict.vertices
meshSpanZ = 0.1
# Diameter of the cylinder
cylDiameter = 0.1
# Reference area
Aref = cylDiameter*meshSpanZ

# Reference density
rhoRef = 1.0

# Dynamic pressure
pDyn = 0.5*rhoRef*(magUInf**2)
forceScaling = 1.0/(Aref*pDyn)

# -------------------------------------------------------------------------- #

set ylabel "$C_d$"
set xrange [0:24]
set yrange [-4.5:4.5]
# $2 corresponds to x component of the fluid force on the cylinder
plot datafile using 1:(forceScaling*$2) with lines title "Drag Coefficient"

set ylabel "$C_l$"
set xrange [0:24]
set yrange [-0.06:0.06]
# $3 corresponds to y component of the fluid force on the cylinder
plot datafile using 1:(forceScaling*$3) with lines title "Lift Coefficient"
