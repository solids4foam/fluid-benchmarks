from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

# Load the uploaded image
img_path = "./Erzincanli2013Fig6.png"
img = Image.open(img_path)

# Convert to grayscale numpy array
gray = img.convert("L")
arr = np.array(gray)

# Invert (so line is bright)
inv = 255 - arr

# Threshold to binary
binary = inv > 100

# Get coordinates of "line" pixels
ys, xs = np.where(binary)

# The plot axes: x from 0 to 10, y from -30 to 30 (approx from the image)
# Map pixel coordinates to data coords
x_min, x_max = 0, 10
y_min, y_max = -30, 30

# Image pixel coords
height, width = arr.shape
# In the image array: y=0 top, so invert y when mapping
x_data = xs / width * (x_max - x_min) + x_min
y_data = (1 - ys / height) * (y_max - y_min) + y_min

# To make cleaner: sort by x and average y in bins
num_bins = 500
bins = np.linspace(x_min, x_max, num_bins)
digitized = np.digitize(x_data, bins)
x_vals = []
y_vals = []
for i in range(1, len(bins)):
    mask = digitized == i
    if np.any(mask):
        x_vals.append(np.mean(x_data[mask]))
        y_vals.append(np.mean(y_data[mask]))

import pandas as pd
df = pd.DataFrame({"Time": x_vals, "Cd": y_vals})

# Export data
df.to_csv("Cd.dat", sep=" ", index=False, header=False)
