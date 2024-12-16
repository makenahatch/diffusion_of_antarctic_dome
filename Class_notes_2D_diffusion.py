#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10

@author: makenahatch

Class notes on 2-D diffusion
"""

### -------- 2-D Diffusion --------

# Solving with FTCS in two dimensions and flattening/reshaping the elevations
# Remember: ik = i * nx + k

import numpy as np
import matplotlib.pyplot as plt

# INITIAL CONDITIONS
nx = 30
ny = 40

dx = 2 # meters
dy = 2 # meters

x = np.arange(0, nx*dx, dx) # creates the 1-D array of x positions
y = np.arange(0, ny*dy, dy) # creates the 1-D array of y positions
X, Y = np.meshgrid(x, y, indexing = 'ij') # creates a 2-D coordinate system for plotting

"""
capital letters = 2D array
lowercase = 1D array
"""

D = 0.02 # m2/year

dt = 5 # years

# example will use topography
Z = np.random.random((nx, ny)) * 100
z = Z.flatten() # gives the 1-D flattened array, flattened by rows

# STABILITY CHECK
sx = dt * D / dx**2
sy = dt * D / dy**2

import sys
if sx > 0.5:
    print('x is unstable')
    sys.exit()
elif sy > 0.5:
    print('y is unstable')
    sys.exit()

# CREATING THE A MATRIX
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

#print(A)

# PLOT INITIAL CONDITIONS
# method 1 - use surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot_surface(X,Y,Z)
ax.set_title('Initial Conditions')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Distance (m)')
ax.set_zlabel('Elevation (m)')

# method 2 = use pcolormesh
fig2, ax2 = plt.subplots(1,1)
initialz = ax2.pcolormesh(X,Y,Z)
ax2.set_title('Initial Conditions')
ax2.set_xlabel('Distance (m)')
ax2.set_ylabel('Distance (m)')
fig2.colorbar(initialz, label = 'Elevation (m)')

# RUN TIME LOOP
time = 0
totaltime = 1000 # years
while time <= totaltime:
    newz = np.dot(A,z)
    z[:] = newz
    time += dt

# PLOT FINAL TOPOGRAPHY
Z = z.reshape(X.shape)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot_surface(X,Y,Z)
ax.set_title('Final Topography')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Distance (m)')
ax.set_zlabel('Elevation (m)')

fig2, ax2 = plt.subplots(1,1)
initialz = ax2.pcolormesh(X,Y,Z)
ax2.set_title('Final Topography')
ax2.set_xlabel('Distance (m)')
ax2.set_ylabel('Distance (m)')
fig2.colorbar(initialz, label = 'Elevation (m)')