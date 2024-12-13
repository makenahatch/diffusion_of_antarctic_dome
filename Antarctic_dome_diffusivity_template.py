#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:53:11 2024

@author: makenahatch
"""

import numpy as np
import matplotlib.pyplot as plt

"""
INITIAL CONDITIONS
"""
nx, ny = 45, 45 # grid size in km/10 (approximates Dome Argus dimensions, dividing by 10 makes the command smaller)
dx, dy = 0.5, 0.5 # grid spacing in km/10

time = 0 # start time in years
totaltime = 1500 # total simulation time in years
dt = 1 # time step in years

# spatial grid
x = np.arange(0, nx*dx, dx) # creates the 1-D array of x positions
y = np.arange(0, ny*dy, dy) # creates the 1-D array of y positions
X, Y = np.meshgrid(x, y, indexing = 'ij') # creates a 2-D coordinate system for plotting

# initial ice thickness in km/10, thicker in the center (approximates Dome Argus profile, dividing by 10 makes the command smaller)
Z_initial = 0.5 * np.exp(-((X - nx * dx / 2) ** 2 + (Y - ny * dy / 2) ** 2) / (2 * (nx * dx / 3) ** 2)) # creates a peak at 0.5 km/10 elevation in the spatial grid 
Z_initial += 0.05  # adds base thickness to simulate broader elevation
z = Z_initial.flatten() # gives the 1-D flattened array, flattened by rows

# diffusivity constant
D = 0.03 # diffusivity parameter for ice flow in (km/10)**2/year (dividing by 10 makes the command smaller)

"""
STABILITY CHECK
"""
sx = dt * D / dx**2
sy = dt * D / dy**2
print(sx,sy)

import sys
if sx > 0.5:
    print('x is unstable')
    sys.exit()
elif sy > 0.5:
    print('y is unstable')
    sys.exit()

"""
CREATING THE A MATRIX
"""
A = np.zeros((nx*ny, nx*ny))

for i in range(nx):
    for k in range(ny):
        ik = i * ny + k
        
        # ---- BOUNDARY CONDITIONS ---- #
        if i == 0:
            A[ik, ik] = 1 # no change
        elif i == (nx-1):
            A[ik, ik] = 1 # no change
        elif k == 0:
            A[ik, ik] = 1 # no change
        elif k == (ny-1):
            A[ik, ik] = 1 # no change
        else:
            # ---- MATRIX COEFFICIENT ---- #
            A[ik, ik] = 1 - 2*sx - 2*sy
            A[ik, (i+1)*ny + k] = sx
            A[ik, (i-1)*ny + k] = sx
            A[ik, i*ny + (k+1)] = sy
            A[ik, i*ny + (k-1)] = sy

"""
RUN TIME LOOP
"""
while time <= totaltime:
    newz = np.dot(A,z)
    z[:] = newz
    time += dt

"""
COMBINED 3-D PLOT
"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot initial conditions
ax.plot_surface(X, Y, Z_initial, color='lightsteelblue', alpha=0.5, label="Initial")

# plot final topography
Z_final = z.reshape(X.shape)
ax.plot_surface(X, Y, Z_final, color='cornflowerblue', alpha=0.8, label="Final")

# labels and title
ax.set_title('Antarctic Ice Sheet Diffusivity Projection')
ax.set_xlabel('Distance (km)')
ax.set_ylabel('Distance (km)')
ax.set_zlabel('Elevation (km)')

# show the plot
plt.legend(['Initial Conditions', 'Final Topography'])
plt.show()