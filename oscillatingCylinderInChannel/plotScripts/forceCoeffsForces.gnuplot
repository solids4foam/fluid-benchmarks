# Plots the drag and the lift coefficients based on the data calculate by the
# `forces` function object.

set terminal pdfcairo enhanced color solid

set datafile separator whitespace
set datafile commentschars "#"
dataFile = "./postProcessing/forces/0/force.dat"
validCdDataFile = "./Cd.dat"
validClDataFile = "./Cl.dat"

set output 'forceCoeffsForces.pdf'

set grid
set tics nomirror
set border 3
set key top left
set key box width 2 height 1

lineColor = "blue"
pointColor = "red"
pointType = 6 # Hollow circle
pointSize = .2

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

# --- Plot Drag Coefficient ---------------------------------------------------#

set title "Drag Coefficient (C_d) vs. Time (t)" . \
          " (using {/Monospace forces} function object)"

set xlabel "t [s]"
set ylabel "C_d [1]"

set xtics 2
set ytics 1.5

set xrange [0:24]
set yrange [-4.5:4.5]

# $2 corresponds to x component of the fluid force on the cylinder
plot \
    dataFile \
        using 1:(forceScaling*$2) \
        title "Present" \
        with lines lc rgb lineColor, \
    validCdDataFile \
        using 1:2 \
        title "Wan\\&Turek" \
        with points pt pointType ps pointSize lc rgb pointColor

# --- Plot Lift Coefficient ---------------------------------------------------#

set title "Lift Coefficient (C_l) vs. Time (t)" . \
          " (using the {/Monospace forces} function object)"

set xlabel "t [s]"
set ylabel "C_l [1]"

set xtics 2
set ytics 0.02

set xrange [0:24]
set yrange [-0.06:0.06]

# $3 corresponds to y component of the fluid force on the cylinder
plot \
    dataFile \
        using 1:(forceScaling*$3) \
        title "Present" \
        with lines lc rgb lineColor, \
    validClDataFile \
        using 1:2 \
        title "Wan\\&Turek" \
        with points pt pointType ps pointSize lc rgb pointColor
