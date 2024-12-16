#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:53:11 2024

@author: makenahatch

This Python script simulates the diffusion of ice sheet thickness of an Antarctic dome, by
approximating the profile of an ice dome in the Antarctic Peninsula region.
"""

# Import libraries
import numpy as np # for numerical operations
import matplotlib.pyplot as plt # for creating plots

### -------- INITIAL CONDITIONS --------

# Form grid
nx, ny = 100, 100 # m
dx, dy = 1, 1 # m
x = np.arange(0, nx*dx, dx) # 1-D array of x positions
y = np.arange(0, ny*dy, dy) # 1-D array of y positions
X, Y = np.meshgrid(x, y, indexing = 'ij') # create 2-D coordinate system for plotting

# Time step
time = 0 # yrs
totaltime = 1500 # yrs
dt = 1 # time step in yrs

# Diffusion
D = 0.2 # m2/yr

# Initial ice thickness
Z_initial = 2000 * np.exp(-((X - nx * dx / 2) ** 2 + (Y - ny * dy / 2) ** 2) / (2 * (nx * dx / 3) ** 2)) # creates a peak at 2,000 m elevation in the spatial grid 
Z_initial += 0.05 # add base thickness to simulate broader elevation
z = Z_initial.flatten() # 1-D flattened array, flattened by rows

# Check stability
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

### -------- CREATE THE A MATRIX --------

A = np.zeros((nx*ny, nx*ny))

for i in range(nx):
    for k in range(ny):
        ik = i * ny + k
        
        # Boundary conditions
        if i == 0:
            A[ik, ik] = 1  # no change
        elif i == (nx-1):
            A[ik, ik] = 1  # no change
        elif k == 0:
            A[ik, ik] = 1  # no change
        elif k == (ny-1):
            A[ik, ik] = 1  # no change
        else:
            # Matrix coefficient
            A[ik, ik] = 1 - 2*sx - 2*sy
            A[ik, (i+1)*ny + k] = sx
            A[ik, (i-1)*ny + k] = sx
            A[ik, i*ny + (k+1)] = sy
            A[ik, i*ny + (k-1)] = sy

### -------- RUN TIME LOOP --------

while time <= totaltime:
    newz = np.dot(A,z)
    z[:] = newz
    time += dt

### --------- PLOT DIFFUSION OF ICE SHEET THICKNESS ---------

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot initial conditions
ax.plot_surface(X, Y, Z_initial, color='lightsteelblue', alpha=0.5, label="Initial ice thickess")

# Plot final topography
Z_final = z.reshape(X.shape)
ax.plot_surface(X, Y, Z_final, color='cornflowerblue', alpha=0.8, label="Final ice thickness")

# Title and labels
ax.set_title('Diffusion of Ice Sheet Thickness')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Distance (m)')
ax.set_zlabel('Elevation (m)')
plt.legend(['Initial ice thickness', 'Final ice thickness'])
plt.show()