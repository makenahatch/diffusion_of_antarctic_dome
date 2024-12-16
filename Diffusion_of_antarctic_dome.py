#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 12:28:06 2024

@author: makenahatch

This script models the diffusion dynamics of an Antarctic dome using surface and bed topography data.
The visualizations generated can be used to help scientists study ice sheet behavior and polar geography.
"""

# Import libraries
import numpy as np # for numerical operations
import matplotlib.pyplot as plt # for creating plots

### -------- INITIAL CONDITIONS --------

# Import .asc data
ascii_grid_surface = np.loadtxt("bedmap_surface.asc", dtype = 'float', skiprows=6)
ascii_grid_bed = np.loadtxt("bedmap_bed.asc", dtype = 'float', skiprows=6)

ascii_headers = np.loadtxt("bedmap_surface.asc", max_rows = 6, dtype = 'str')
n_long = ascii_headers[0,1].astype(int)
n_lat = ascii_headers[1,1].astype(int)
dxy = ascii_headers[4,1].astype(float)
xllcorner = ascii_headers[2,1].astype(float)
yllcorner = ascii_headers[3,1].astype(float)

# Form grid
x = np.arange(0, dxy*n_lat, dxy) + xllcorner # 1-D array of x positions
y = np.arange(0, dxy*n_long, dxy) + yllcorner # 1-D array of y positions
LAT, LONG = np.meshgrid(x, y, indexing='ij') # create 2-D coordinate system for plotting
nodes = n_long*n_lat

# Time step
time = 0 # yrs
totaltime = 100000 # yrs
dt = 100 # time step in yrs

"""
-------- ARTIFICIAL DIFFUSION --------

This process is driven by advective properties, but makes something that we intuitively want to model
as a diffusive process, so we call it artificial diffusion.

Publication introducing artificial diffusion:
    MacAyeal, D. R. EISMINT: Lessons in Ice-Sheet Modeling. Department of Geophysical Sciences,
    University of Chicago.

Ice velocity data:
    Value used for ice velocity below is estimated using the linked maps. Value may be slightly
    underestimated due to changing the location of interest.
    https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-ice-velocity-maps
    https://nsidc.org/aiv
"""

# Artificial diffusion
ice_sheet_velocity = 2 # m/yr
D = ice_sheet_velocity * dxy / 2 # artificial diffusion rate in m2/yr

# Check stability
sx = dt * D / dxy**2
sy = dt * D / dxy**2
print(sx,sy) 

import sys
if sx > 0.5:
    print('x is unstable')
    sys.exit()
elif sy > 0.5:
    print('y is unstable')
    sys.exit()

### -------- CREATE THE A MATRIX --------

A = np.zeros((n_long*n_lat, n_long*n_lat))

for i in range(n_lat):
    for k in range(n_long):
        ik = i * n_long + k
        
        # Boundary conditions
        if i == 0:
            A[ik, ik] = 1  # no change
        elif i == (n_lat-1):
            A[ik, ik] = 1  # no change
        elif k == 0:
            A[ik, ik] = 1  # no change
        elif k == (n_long-1):
            A[ik, ik] = 1  # no change
        else:
            # Matrix coefficient
            A[ik, ik] = 1 - 2*sx - 2*sy
            A[ik, (i+1)*n_long + k] = sx
            A[ik, (i-1)*n_long + k] = sx
            A[ik, i*n_long + (k+1)] = sy
            A[ik, i*n_long + (k-1)] = sy

### -------- RUN TIME LOOP --------

# Initial elevation
elv_flat = ascii_grid_surface.flatten() # 1-D flattened array, flattened by rows

while time <= totaltime:
    new_elv = np.dot(A,elv_flat)
    elv_flat[:] = new_elv
    time += dt
   
### --------- PLOT DIFFUSION OF ICE SHEET ---------

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot initial surface elevation
ax.plot_surface(LAT, LONG, ascii_grid_surface, color='royalblue', alpha=0.3)

# Plot surface elevation after 100,000 yrs
final_elv = elv_flat.reshape(LAT.shape)
ax.plot_surface(LAT, LONG, final_elv, color='cornflowerblue', alpha=0.8)

# Plot bedrock elevation
ax.plot_surface(LAT, LONG, ascii_grid_bed, color='slategrey', alpha=0.3)

# Title and labels
ax.set_title('100,000-Year Diffusion of Ice Sheet') 
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Distance (m)')
ax.set_zlabel('Elevation (m)')
plt.legend(['Initial surface elevation', '100,000-year surface elevation', 'Bedrock elevation'])
plt.show()

### --------- PLOT CHANGE IN ICE THICKNESS ---------

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot change in ice thickness
initial_ice_thickness = ascii_grid_surface - ascii_grid_bed
final_ice_thickness = final_elv - ascii_grid_bed
change_in_ice_thickness = initial_ice_thickness - final_ice_thickness
ax.plot_surface(LAT, LONG, change_in_ice_thickness, color='mediumpurple', alpha=0.3)

# Title and labels
ax.set_title('100,000-Year Change in Ice Thickness')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Distance (m)')
ax.set_zlabel('Distance (m)')
plt.show()

### --------- SAVE ASCII OUTPUT ---------

header = 'NCOLS %s \n' % n_long + 'NROWS %s \n' % n_lat + 'xllcorner %s \n' % xllcorner+ 'yllcorner %s \n' % yllcorner + 'cellsize %s \n' % dxy + 'NODATA_value -9999'
np.savetxt('new_elev.asc', final_elv, header = header, comments = '')