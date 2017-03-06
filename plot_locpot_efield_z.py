# -*- coding: utf-8 -*-
#!/usr/bin/python

import pyximport; pyximport.install()
import sys
import numpy as np
import locpot_import
import periodic_grid_translation
import periodic_atom_translation
from nonorthogonal_to_orthogonal_grid_interpolation import interpolate_using_linear_transformation, interp_linear_transform_orthogonal_z
from plot_dft_grid_data import plot_efield, plot_pot_yline, plot_pot_zline, plot_Ey_yline, plot_Ez_zline
from plot_atoms import plot_atoms
import matplotlib.pyplot as plt


def plot_efield_z(file_path, grid_spacing, plot_pos, plot_box_size, efield_xy_limits, efield_xz_limits, atom_radius_scaling):
    is_grid_orthogonal, is_z_orthogonal = locpot_import.check_grid_orthogonality(file_path)
    
    if is_grid_orthogonal:
        print 'Reading LOCPOT file with orthogonal grid'
        atom_types, atom_positions, xs, ys, zs, pot_on_grid = locpot_import.read_locpot(file_path, orthogonal=True)
        print 'nx = {}, ny = {}, nz = {}'.format(xs.shape[0], ys.shape[0], zs.shape[0])
        anchor_atom = np.argmax(atom_positions[:, 2])
        anchor_atom_pos = atom_positions[anchor_atom, :]
        print 'Translating the grid'
        xs, ys, zs, pot_on_grid = periodic_grid_translation.translate_grid(xs, ys, zs, pot_on_grid, anchor_atom_pos)
    else:
        print 'Reading LOCPOT file non-orthogonal grid'
        atom_types, atom_positions, n_voxels, voxel_vectors, nonorthogonal_pot_grid = \
            locpot_import.read_locpot(file_path, orthogonal=False, return_unstructured_data=False)
        print 'nx_voxels = {}, ny_voxels = {}, nz_voxels = {}'.format(n_voxels[0], n_voxels[1], n_voxels[2])
        anchor_atom = np.argmax(atom_positions[:, 2])
        anchor_atom_pos = atom_positions[anchor_atom, :]
        print 'Interpolating to regular, orthogonal grid'
        if is_z_orthogonal:
            xs, ys, zs, pot_on_grid = interp_linear_transform_orthogonal_z(nonorthogonal_pot_grid,
                                            n_voxels, voxel_vectors, anchor_atom_pos, grid_spacing)
        else:
            xs, ys, zs, pot_on_grid = interpolate_using_linear_transformation(nonorthogonal_pot_grid,
                                            n_voxels, voxel_vectors, anchor_atom_pos, grid_spacing)
        print 'nx = {}, ny = {}, nz = {}'.format(xs.shape[0], ys.shape[0], zs.shape[0])
    pot_on_grid = -pot_on_grid

    # Translate the coordinate system so that the tip atom is at position (0, 0, 0)
    xs = xs - anchor_atom_pos[0]
    ys = ys - anchor_atom_pos[1]
    zs = zs - anchor_atom_pos[2]
    if is_grid_orthogonal:
        n_voxels = np.array([xs.shape[0], ys.shape[0], zs.shape[0]])
        voxel_vectors = np.array([[xs[1]-xs[0], 0,           0],
                                  [0,           ys[1]-ys[0], 0],
                                  [0,           0,           zs[1]-zs[0]]])
    atom_positions = periodic_atom_translation.translate_atoms(atom_positions, n_voxels, voxel_vectors, -anchor_atom_pos)
    
    print 'Plotting'
    
    ax1 = plt.subplot(2, 1, 1)
    plot_efield(xs, ys, zs, pot_on_grid, plot_pos[1], slice_normal='y', efield_component='z', limits=efield_xz_limits)
    plot_atoms(atom_types, atom_positions, normal_direction='y', atom_radius_scaling=atom_radius_scaling)
    plt.xlim(-plot_box_size/2, plot_box_size/2)
    plt.ylim(-plot_box_size/3, 2*plot_box_size/3)
    ax1.set_aspect('equal')
    plt.title(u'Electric field (cut through tip)', size=16)
    plt.xlabel(u'x (Å)', size=14)
    plt.ylabel(u'z (Å)', size=14)

    ax2 = plt.subplot(2, 1, 2)
    plot_efield(xs, ys, zs, pot_on_grid, plot_pos[2], slice_normal='z', efield_component='z', limits=efield_xy_limits)
    plot_atoms(atom_types, atom_positions, normal_direction='z', height_cut=0.0, atom_radius_scaling=atom_radius_scaling, atom_alpha=0.4)
    plt.xlim(-plot_box_size/2, plot_box_size/2)
    plt.ylim(-plot_box_size/2, plot_box_size/2)
    ax2.set_aspect('equal')
    plt.title(u'Electric field ({} Å above surface)'.format(plot_pos[2]), size=16)
    plt.xlabel(u'x (Å)', size=14)
    plt.ylabel(u'y (Å)', size=14)
