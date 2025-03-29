# Copyright (c) 2025 Vito Romanelli Tricanico
#
# SPDX-License-Identifier: BSD-2-Clause (see LICENSE)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

data = np.loadtxt("output.txt")

# Create figure and axes handler with black background
plt.style.use('dark_background')
plt.rcParams['grid.color'] = "black" # grid lines are effectivelly black (for 3d plots)
fig = plt.figure()
fig.patch.set_facecolor('black')  # Set figure background color to black
ax = fig.add_subplot(111, projection="3d")

# Make panes transparent 
ax.xaxis.pane.fill = False # x pane
ax.yaxis.pane.fill = False # y pane
ax.zaxis.pane.fill = False # z pane

# Set edges to black

ax.xaxis.pane.set_edgecolor('black')
ax.yaxis.pane.set_edgecolor('black')
ax.zaxis.pane.set_edgecolor('black')
ax.xaxis.line.set_visible(False)
ax.yaxis.line.set_visible(False)
ax.zaxis.line.set_visible(False)

# Set ticks (numbers) to black

#ax.tick_params(axis='x', colors='black')
#ax.tick_params(axis='y', colors='black')
#ax.tick_params(axis='z', colors='black')

n_values = data[:, 0]
x_values = data[:, 1]
y_values = data[:, 2]
z_values = data[:, 3]

#print(np.unique(n_values))

# animation
def update(frame):
    ax.cla() # Clear current axes (comment if plotting paths is desired)
    # set axes limits
    
    #ax.set_xlim3d(-20, 40)
    #ax.set_ylim3d(-20, 40)
    #ax.set_zlim3d(-20, 40)
    
    idx = np.where(n_values == frame)
    return  ax.scatter(x_values[idx], y_values[idx], z_values[idx], c="greenyellow", s=5, edgecolor="black", linewidth=0)

# Create the animation
ani = FuncAnimation(fig, update, frames = np.unique(n_values))
# Save the animation as an MP4 file
ani.save("animation.mp4", writer="ffmpeg", fps=20)

plt.close(fig)
