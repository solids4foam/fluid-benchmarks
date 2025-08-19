#!/usr/bin/env python3

import numpy as np

# Parameters
A = 0.25        # amplitude
f = 0.25        # frequency [Hz]
t0 = 0.0        # start time
end_time = 25.0 # end time
delta_t = 0.005 # time step
x0 = 0.0        # initial displacement

# Time vector
time = np.arange(t0, end_time + delta_t, delta_t)

# Displacement (only x-component varies)
displacement = x0 + A * np.sin(2 * np.pi * f * time)

# Write output
outfile = "timeVsDisplacement"
with open(outfile, "w") as file:
    file.write("(\n")
    for ti, xi in zip(time, displacement):
        file.write(f"    ( {ti:.3f} ( {xi:.6f} 0 0 ))\n")
    file.write(")\n")

print(f"Data written to '{outfile}' with {len(time)} entries.")
