# -*- coding: utf-8 -*-
#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

divergent_colormap = plt.cm.bwr
sequential_colormap = plt.cm.gnuplot2
n_contour_default = 80

def interpolate_slice_from_grid(xs, ys, zs, data_grid, pos_along_normal, normal_direction):
    if normal_direction == 'x':
        normal_coord = 0
        normal_grid_spacing = xs[1] - xs[0]
        normal_origin = xs[0]
    elif normal_direction == 'y':
        normal_coord = 1
        normal_grid_spacing = ys[1] - ys[0]
        normal_origin = ys[0]
    else:
        normal_coord = 2
        normal_grid_spacing = zs[1] - zs[0]
        normal_origin = zs[0]
    
    i_normal_lower = int((pos_along_normal - normal_origin)/normal_grid_spacing)
    i_normal_higher = i_normal_lower + 1
    interp_factor = (pos_along_normal - normal_origin)/normal_grid_spacing - i_normal_lower
    
    if normal_coord == 0:
        data_slice = np.zeros((ys.shape[0], zs.shape[0]))
        for iy in range(data_slice.shape[0]):
            for iz in range(data_slice.shape[1]):
                data_slice[iy, iz] = (1-interp_factor)*data_grid[i_normal_lower, iy, iz] + interp_factor*data_grid[i_normal_higher, iy, iz]
        return ys, zs, data_slice
    
    elif normal_coord == 1:
        data_slice = np.zeros((xs.shape[0], zs.shape[0]))
        for ix in range(data_slice.shape[0]):
            for iz in range(data_slice.shape[1]):
                data_slice[ix, iz] = (1-interp_factor)*data_grid[ix, i_normal_lower, iz] + interp_factor*data_grid[ix, i_normal_higher, iz]
        return xs, zs, data_slice
    
    else:
        data_slice = np.zeros((xs.shape[0], ys.shape[0]))
        for ix in range(data_slice.shape[0]):
            for iy in range(data_slice.shape[1]):
                data_slice[ix, iy] = (1-interp_factor)*data_grid[ix, iy, i_normal_lower] + interp_factor*data_grid[ix, iy, i_normal_higher]
        return xs, ys, data_slice


def plot_density(xs, ys, zs, density, slice_pos, slice_normal, limits=None, n_contour=n_contour_default, 
                contour_levels=np.array([]), add_colorbar=True, colormap='sequential'):
    if colormap=='divergent':
        cmap = divergent_colormap
    else:
        cmap = sequential_colormap
    
    #interpolate slice
    coords_1, coords_2, dens_slice = interpolate_slice_from_grid(xs, ys, zs, density, slice_pos, normal_direction=slice_normal)
    
    #plotting
    if contour_levels.size == 0 and limits is None:
        cnt = plt.contourf(coords_1, coords_2, dens_slice.T, n_contour, cmap=cmap, alpha=1.0)
    else:
        if contour_levels.size == 0:
            dens_slice = dens_slice.clip(limits[0], limits[1])
            contour_levels = np.linspace(limits[0], limits[1], n_contour, endpoint=True)
        cnt = plt.contourf(coords_1, coords_2, dens_slice.T, contour_levels, cmap=cmap, alpha=1.0)
    
    # This is the fix for the white lines between contour levels
    for c in cnt.collections:
        c.set_edgecolor("face")
    
    if add_colorbar:
        cbar = plt.colorbar() #plt.colorbar(pad=0.05, fraction=0.05, shrink=0.9)
        cbar.ax.set_ylabel(u'Charge density (e/Å$^3$)')
        if limits is not None:
            cbar.set_ticks(np.linspace(limits[0], limits[1], 5, endpoint=True))


def plot_pot(xs, ys, zs, pot_on_grid, slice_pos, slice_normal, limits=None, n_contour=n_contour_default,
            contour_levels=np.array([]), add_colorbar=True, colormap='divergent'):
    if colormap=='divergent':
        cmap = divergent_colormap
    else:
        cmap = sequential_colormap
    
    #interpolate slice
    coords_1, coords_2, pot_slice = interpolate_slice_from_grid(xs, ys, zs, pot_on_grid, slice_pos, normal_direction=slice_normal)
    
    #plotting
    if contour_levels.size == 0 and limits is None:
        cnt = plt.contourf(coords_1, coords_2, pot_slice.T, n_contour, cmap=cmap, alpha=1.0)
    else:
        if contour_levels.size == 0:
            pot_slice = pot_slice.clip(limits[0], limits[1])
            contour_levels = np.linspace(limits[0], limits[1], n_contour, endpoint=True)
        cnt = plt.contourf(coords_1, coords_2, pot_slice.T, contour_levels, cmap=cmap, alpha=1.0)
    
    # This is the fix for the white lines between contour levels
    for c in cnt.collections:
        c.set_edgecolor("face")
    
    if add_colorbar:
        cbar = plt.colorbar() #plt.colorbar(pad=0.05, fraction=0.05, shrink=0.9)
        cbar.ax.set_ylabel(u'$\Phi$ (V)')
        if limits is not None:
            cbar.set_ticks(np.linspace(limits[0], limits[1], 5, endpoint=True))


def plot_efield(xs, ys, zs, pot_on_grid, slice_pos, slice_normal, efield_component, limits=None, n_contour=n_contour_default,
                contour_levels=np.array([]), add_colorbar=True, colormap='divergent'):
    if colormap=='divergent':
        cmap = divergent_colormap
    else:
        cmap = sequential_colormap
    
    #calculate the electric field
    dx = xs[1]-xs[0]
    dy = ys[1]-ys[0]
    dz = zs[1]-zs[0]
    Ex, Ey, Ez = np.gradient(-pot_on_grid, dx, dy, dz)
    
    if efield_component == 'x':
        efield = Ex
    elif efield_component == 'y':
        efield = Ey
    else:
        efield = Ez
    
    #interpolate slice
    coords_1, coords_2, efield_slice = interpolate_slice_from_grid(xs, ys, zs, efield, slice_pos, normal_direction=slice_normal)
    
    #plotting
    if contour_levels.size == 0 and limits is None:
        cnt = plt.contourf(coords_1, coords_2, efield_slice.T, n_contour, cmap=cmap, alpha=1.0)
    else:
        if contour_levels.size == 0:
            efield_slice = efield_slice.clip(limits[0], limits[1])
            contour_levels = np.linspace(limits[0], limits[1], n_contour, endpoint=True)
        cnt = plt.contourf(coords_1, coords_2, efield_slice.T, contour_levels, cmap=cmap, alpha=1.0)
    
    # This is the fix for the white lines between contour levels
    for c in cnt.collections:
        c.set_edgecolor("face")
    
    if add_colorbar:
        cbar = plt.colorbar() #plt.colorbar(pad=0.05, fraction=0.05, shrink=0.9)
        cbar.ax.set_ylabel(u'$E_{}$ (V/Å)'.format(efield_component))
        if limits is not None:
            cbar.set_ticks(np.linspace(limits[0], limits[1], 5, endpoint=True))


def plot_dzEz_xy(xs, ys, zs, pot_on_grid, z_slice):
    #calculate the gradient of the z component of electric field
    dx = xs[1]-xs[0]
    dy = ys[1]-ys[0]
    dz = zs[1]-zs[0]
    Ex, Ey, Ez = np.gradient(-pot_on_grid, dx, dy, dz)
    dxEz, dyEz, dzEz = np.gradient(Ez, dx, dy, dz)
    
    #interpolate slice
    xs, ys, dzEz_slice = interpolate_slice_from_grid(xs, ys, zs, dzEz, z_slice, normal_direction='z')
    
    #plotting
    plt.contourf(xs, ys, dzEz_slice.T, n_contour, cmap=colormap, alpha=1.0)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(u'$dE_z/dz$ (V/Å$^2$)', size=14)


def plot_pot_yline(xs, ys, zs, pot_on_grid, xz_point, y_min, y_max, n_points=100):
    pot_interpolator = RegularGridInterpolator((xs, ys, zs), pot_on_grid)
    
    #interpolate pot to line
    ys_line = np.linspace(y_min, y_max, num=n_points)
    yline_points = np.ones((len(ys_line), 3))
    yline_points[:, 0] = yline_points[:, 0]*xz_point[0]
    yline_points[:, 1] = ys_line
    yline_points[:, 2] = yline_points[:, 2]*xz_point[1]
    pot_line = pot_interpolator(yline_points)
    
    #plotting
    plt.plot(ys_line, pot_line)
    
    return ys_line, pot_line


def plot_pot_zline(xs, ys, zs, pot_on_grid, xy_point, z_min, z_max, n_points=100):
    pot_interpolator = RegularGridInterpolator((xs, ys, zs), pot_on_grid)
    
    #interpolate pot to line
    zs_line = np.linspace(z_min, z_max, num=n_points)
    zline_points = np.ones((len(zs_line), 3))
    zline_points[:, 0] = zline_points[:, 0]*xy_point[0]
    zline_points[:, 1] = zline_points[:, 1]*xy_point[1]
    zline_points[:, 2] = zs_line
    pot_line = pot_interpolator(zline_points)
    
    #plotting
    plt.plot(zs_line, pot_line)
    
    return zs_line, pot_line


def plot_Ey_yline(xs, ys, zs, pot_on_grid, xz_point, y_min, y_max, n_points=100):
    #calculate the electric field
    dx = xs[1]-xs[0]
    dy = ys[1]-ys[0]
    dz = zs[1]-zs[0]
    Ex, Ey, Ez = np.gradient(-pot_on_grid, dx, dy, dz, edge_order=2)
    
    Ey_interpolator = RegularGridInterpolator((xs, ys, zs), Ey)
    
    #interpolate Ey to line
    ys_line = np.linspace(y_min, y_max, num=n_points)
    yline_points = np.ones((len(ys_line), 3))
    yline_points[:, 0] = yline_points[:, 0]*xz_point[0]
    yline_points[:, 1] = ys_line
    yline_points[:, 2] = yline_points[:, 2]*xz_point[1]
    Ey_line = Ey_interpolator(yline_points)
    
    #plotting
    plt.plot(ys_line, Ey_line)


def plot_Ez_zline(xs, ys, zs, pot_on_grid, xy_point, z_min, z_max, n_points=100):
    #calculate the electric field
    dx = xs[1]-xs[0]
    dy = ys[1]-ys[0]
    dz = zs[1]-zs[0]
    Ex, Ey, Ez = np.gradient(-pot_on_grid, dx, dy, dz)
    
    Ez_interpolator = RegularGridInterpolator((xs, ys, zs), Ez)
    
    #interpolate Ez to line
    zs_line = np.linspace(z_min, z_max, num=n_points)
    zline_points = np.ones((len(zs_line), 3))
    zline_points[:, 0] = zline_points[:, 0]*xy_point[0]
    zline_points[:, 1] = zline_points[:, 1]*xy_point[1]
    zline_points[:, 2] = zs_line
    Ez_line = Ez_interpolator(zline_points)
    
    #plotting
    plt.plot(zs_line, Ez_line)


def plot_dzEz_zline(xs, ys, zs, pot_on_grid, xy_point, z_min, z_max, n_points=100):
    #calculate the electric field
    dx = xs[1]-xs[0]
    dy = ys[1]-ys[0]
    dz = zs[1]-zs[0]
    Ex, Ey, Ez = np.gradient(-pot_on_grid, dx, dy, dz)
    dxEz, dyEz, dzEz = np.gradient(Ez, dx, dy, dz)
    
    dzEz_interpolator = RegularGridInterpolator((xs, ys, zs), dzEz)
    
    #interpolate Ez to line
    zs_line = np.linspace(z_min, z_max, num=n_points)
    zline_points = np.ones((len(zs_line), 3))
    zline_points[:, 0] = zline_points[:, 0]*xy_point[0]
    zline_points[:, 1] = zline_points[:, 1]*xy_point[1]
    zline_points[:, 2] = zs_line
    dzEz_line = dzEz_interpolator(zline_points)
    
    #plotting
    plt.plot(zs_line, dzEz_line)
