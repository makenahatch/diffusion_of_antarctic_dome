#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 12:28:06 2024

@author: makenahatch

This Python script models and visualizes diffusivity over a region near Dome Argus in Antarctica,
based on BEDMAP2 data. It begins by importing necessary libraries (numpy for numerical operations
and matplotlib for plotting) and reading input data from an ASCII file. The script processes this
file to extract geographic grid information (such as dimensions, resolution, and corner
coordinates) and sets up a 2D grid for spatial calculations. Initial conditions, including time
variables and step size, are defined for a simulation of diffusion processes over the spatial
grid, representing phenomena like changes in ice thickness.
"""

import numpy as np
import matplotlib.pyplot as plt

### -------- INITIAL CONDITIONS -------------

# import .asc data
ascii_grid = np.loadtxt("dome_A_thickness.asc", dtype = 'float', skiprows=6)

ascii_headers = np.loadtxt("dome_A_thickness.asc", max_rows = 6, dtype = 'str')
n_long = ascii_headers[0,1].astype(int)
n_lat = ascii_headers[1,1].astype(int)
dxy = ascii_headers[4,1].astype(float)
xllcorner = ascii_headers[2,1].astype(float)
yllcorner = ascii_headers[3,1].astype(float)

## ---- GRID FORMATION
x = np.arange(0, dxy*n_lat, dxy) + xllcorner # array of x values
y = np.arange(0, dxy*n_long, dxy) + yllcorner # array of z values
LAT, LONG = np.meshgrid(x, y, indexing='ij') # this sets up a plotting grid
nodes = n_long*n_lat

## ---- TIME STEP
time = 0 # start time in years
totaltime = 10 # total simulation time in years
dt = 0.1 # time step in years

## ---- DIFFUSION
D = 2600000 # m2/year

## ---- INITIAL ELEVATION
elv_flat = ascii_grid.flatten()

## ---- STABILITY CHECK
sx = dt * D / dxy**2
sy = dt * D / dxy**2
print(sx,sy)

import sys
if sx > 0.5:
    print('n_lat is unstable')
    sys.exit()
elif sy > 0.5:
    print('n_long is unstable')
    sys.exit()

### ---------- A MATRIX AND RUNNING THE MODEL -------------

A = np.zeros((n_long*n_lat, n_long*n_lat))

for i in range(n_lat):
    for k in range(n_long):
        ik = i * n_long + k
        
        # ---- BOUNDARY CONDITIONS ----
        if i == 0:
            A[ik, ik] = 1 # no change
        elif i == (n_lat-1):
            A[ik, ik] = 1 # no change
        elif k == 0:
            A[ik, ik] = 1 # no change
        elif k == (n_long-1):
            A[ik, ik] = 1 # no change
        else:
            # ---- MATRIX COEFFICIENT ----
            A[ik, ik] = 1 - 2*sx - 2*sy
            A[ik, (i+1)*n_long + k] = sx
            A[ik, (i-1)*n_long + k] = sx
            A[ik, i*n_long + (k+1)] = sy
            A[ik, i*n_long + (k-1)] = sy

while time <= totaltime:
    new_elv = np.dot(A,elv_flat)
    elv_flat[:] = new_elv
    time += dt
   
### --------- COMBINED 3-D PLOT ---------

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot initial conditions
ax.plot_surface(LAT, LONG, ascii_grid, color='lightsteelblue', alpha=0.5, label="Initial")

# plot final topography
final_elv = elv_flat.reshape(LAT.shape)
ax.plot_surface(LAT, LONG, final_elv, color='cornflowerblue', alpha=0.8, label="Final")

# labels and title
ax.set_title('Dome Argus Diffusivity Projection')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Distance (m)')
ax.set_zlabel('Elevation (m)')

# show the plot
plt.legend(['Initial Conditions', 'Final Topography'])
plt.show()

### --------- SAVE ASCII OUTPUT ---------

header = 'NCOLS %s \n' % n_long + 'NROWS %s \n' % n_lat + 'xllcorner %s \n' % xllcorner+ 'yllcorner %s \n' % yllcorner + 'cellsize %s \n' % dxy + 'NODATA_value -9999'
np.savetxt('new_elev.asc', final_elv, header = header, comments = '')