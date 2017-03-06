# -*- coding: utf-8 -*-
#!/usr/bin/python

import numpy as np
import scipy.interpolate as sp_interpolate

# Works for data without grid structure, i.e. pot_data having rows
# x y z value
# Interpolates the data into a regular grid.

def interpolate_to_grid_nonortho_z(pot_data, grid_spacing):
    #create suitable orthogonal grid
    xmin = np.amin(pot_data[:, 0])
    xmax = np.amax(pot_data[:, 0])
    nx = int((xmax-xmin)/grid_spacing)
    xs = np.linspace(xmin, xmax, nx)
    
    ymin = np.amin(pot_data[:, 1])
    ymax = np.amax(pot_data[:, 1])
    ny = int((ymax-ymin)/grid_spacing)
    ys = np.linspace(ymin, ymax, ny)
    
    zmin = np.amin(pot_data[:, 2])
    zmax = np.amax(pot_data[:, 2])
    nz = int((zmax-zmin)/grid_spacing)
    zs = np.linspace(zmin, zmax, nz)
    
    xgrid, ygrid, zgrid = np.meshgrid(xs, ys, zs)
    
    #interpolate electrostatic potential to the orthogonal grid
    pot_on_grid = sp_interpolate.griddata((pot_data[:, 0], pot_data[:, 1], pot_data[:, 2]), pot_data[:, 3], (xgrid, ygrid, zgrid), method='linear')
    
    return xs, ys, zs, pot_on_grid


def interpolate_to_grid_orthogonal_z(pot_data, grid_spacing, nz):
    #create suitable orthogonal grid
    xmin = np.amin(pot_data[:, 0])
    xmax = np.amax(pot_data[:, 0])
    nx = int((xmax-xmin)/grid_spacing)
    xs = np.linspace(xmin, xmax, nx)
    
    ymin = np.amin(pot_data[:, 1])
    ymax = np.amax(pot_data[:, 1])
    ny = int((ymax-ymin)/grid_spacing)
    ys = np.linspace(ymin, ymax, ny)
    
    zs = pot_data[0:nz, 2]
    
    xgrid, ygrid = np.meshgrid(xs, ys, indexing='ij')
    
    pot_on_grid = np.zeros((nx, ny, nz))
    #for each z interpolate electrostatic potential to the orthogonal grid in xy
    for iz in range(nz):
        if (iz+1)%10 == 0:
            print 'z layer ' + repr(iz+1) + '/' + repr(nz)
        pot_xy = pot_data[iz::nz, 3]
        pot_on_grid[:, :, iz] = sp_interpolate.griddata((pot_data[iz::nz, 0], pot_data[iz::nz, 1]),
                                                        pot_xy, (xgrid, ygrid), method='linear')
    return xs, ys, zs, pot_on_grid
    
