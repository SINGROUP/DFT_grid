# -*- coding: utf-8 -*-

import pyximport; pyximport.install()
import numpy as np
import cube_import
from nonorthogonal_to_orthogonal_grid_interpolation import interpolate_using_linear_transformation
from plot_dft_grid_data import plot_pot, plot_pot_yline, plot_pot_zline, plot_Ey_yline, plot_Ez_zline
from plot_atoms import plot_atoms
import matplotlib.pyplot as plt

atom_radius_scaling = 0.3

def plot_epot_y_cube(cube_file, description, fig_format, is_nonorthorhombic, efield_relative_pos, grid_spacing=0.1,
                    pot_limits=None, efield_limits=None, do_plot_atoms=False, rel_atom_plot_limits=[0.0, 1.0], colormap='sequential'):

    if is_nonorthorhombic:
        print 'Reading cube file with non-orthogonal grid: {}'.format(cube_file)
        atom_types, atom_pos, n_voxels, voxel_vectors, pot_on_nonortho_grid = cube_import.read_cube_nonorthogonal(cube_file, return_unstructured_data = False)
        print 'nx_voxels = {}, ny_voxels = {}, nz_voxels = {}'.format(n_voxels[0], n_voxels[1], n_voxels[2])
        print 'Interpolating to regular, orthogonal grid'
        new_origin = np.array([0.0, 0.0, 0.0])
        xs, ys, zs, pot_on_grid = interpolate_using_linear_transformation(pot_on_nonortho_grid,
                                                n_voxels, voxel_vectors, new_origin, grid_spacing)
        nx = xs.shape[0]
        ny = ys.shape[0]
        nz = zs.shape[0]
    else:
        print 'Reading cube file with orthogonal grid: {}'.format(cube_file)
        atom_types, atom_pos, xs, ys, zs, pot_on_grid = cube_import.read_cube_orthogonal(cube_file)
        nx = len(xs)
        ny = len(ys)
        nz = len(zs)
    
    print 'nx = {}, ny = {}, nz = {}'.format(nx, ny, nz)

    # Definition of Hartree potential
    pot_on_grid = -pot_on_grid

    print 'Plotting'

    cell_x = (xs[1]-xs[0]) * nx
    cell_y = (ys[1]-ys[0]) * ny
    cell_z = (zs[1]-zs[0]) * nz
    efield_x = efield_relative_pos * cell_x + xs[0]
    efield_y = efield_relative_pos * cell_y + ys[0]
    efield_z = efield_relative_pos * cell_z + zs[0]
    atom_plot_limits = np.array(rel_atom_plot_limits) * cell_z + zs[0]

    plt.figure(figsize=(11, 10))
    plot_pot(xs, ys, zs, pot_on_grid, efield_z, slice_normal='z', limits=pot_limits, colormap=colormap)
    if do_plot_atoms:
        plot_atoms(atom_types, atom_pos, normal_direction='z', atom_plot_limits=atom_plot_limits, atom_radius_scaling=atom_radius_scaling)
    plt.axes().set_aspect('equal')
    plt.title(u'Electrostatic potential')
    plt.xlabel(u'x (Å)')
    plt.ylabel(u'y (Å)')
    plt.tight_layout()
    print 'Saving figure to file {}'.format(description + '-pot_xy' + fig_format)
    plt.savefig(description + '-pot_xy' + fig_format, bbox_inches='tight')

    plt.figure(figsize=(10, 15))
    plt.subplot(211)
    ys_line, pot_line = plot_pot_yline(xs, ys, zs, pot_on_grid, np.array([efield_x, efield_z]), ys[0], ys[-1], n_points=ny)
    plt.xlabel(u'y (Å)')
    plt.ylabel(u'$\Phi$ (V)')
    plt.ylim(pot_limits[0], pot_limits[1])
    plt.subplot(212)
    plot_Ey_yline(xs, ys, zs, pot_on_grid, np.array([efield_x, efield_z]), ys[0], ys[-1], n_points=ny)
    plt.xlabel(u'y (Å)')
    plt.ylabel(u'$E_y$ (V/Å)')
    plt.ylim(efield_limits[0], efield_limits[1])
    plt.tight_layout()
    print 'Saving figure to file {}'.format(description + '-pot_efield_y' + fig_format)
    plt.savefig(description + '-pot_efield_y' + fig_format, bbox_inches='tight')
    
    print 'Done'
