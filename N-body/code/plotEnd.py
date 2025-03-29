# Copyright (c) 2025 Vito Romanelli Tricanico
#
# SPDX-License-Identifier: BSD-2-Clause (see LICENSE)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter

# Load data
data = np.loadtxt("output.txt")

# Create figure and axes handler with black background
plt.style.use('dark_background')
fig = plt.figure()
fig.patch.set_facecolor('black')  # Set figure background color to black
ax = fig.add_subplot(111, projection="3d")

# Remove all axes
ax.set_axis_off()

# Extract data
n_values = data[:, 0]
x_values = data[:, 1]
y_values = data[:, 2]
z_values = data[:, 3]

# Scatter plot of the last frame of data
uni = np.unique(n_values)
idx = np.where(n_values == uni[-1])
ax.scatter(x_values[idx], y_values[idx], z_values[idx], c="blueviolet", s=7, edgecolor="black", linewidth=0)

# Define camera rotation function
def update(frame):
    ax.view_init(elev=30, azim=frame)  # Change elevation and azimuth angles

# Set up animation for rotating the camera
ani = FuncAnimation(fig, update, frames=np.arange(0, 360, 2), interval=50)

# Save the animation as a GIF file using PillowWriter
output_file = "end.gif"
writer = PillowWriter(fps=20)
ani.save(output_file, writer=writer)

print(f"Animation saved as {output_file}")

