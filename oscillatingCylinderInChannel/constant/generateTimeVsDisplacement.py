import numpy as np

# Parameters
A = 0.25
f = 0.25
dt = 0.005
endtime = 24

# Time array from 0 to endtime inclusive
t_vals = np.arange(0, 1 + endtime, dt)

# Compute z(t) = A * sin(2πft)
z_vals = A * np.sin(2 * np.pi * f * t_vals)

# Format the output as required
formatted_list = ["( {:.3f} ( {:.6f} 0 0 ))".format(t, z) for t, z in zip(t_vals, z_vals)]

# Add outer parentheses and join with newlines
output_string = "(\n" + "\n".join(formatted_list) + "\n)"

# Write to file
with open("timeVsDisplacement", "w") as file:
    file.write(output_string)

print("File 'timeVsDisplacement' has been written.")
