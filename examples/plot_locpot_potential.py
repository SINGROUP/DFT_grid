# -*- coding: utf-8 -*-
#!/usr/bin/python

import pyximport; pyximport.install()
import sys
import numpy as np
import matplotlib.pyplot as plt

import locpot_import
import periodic_grid_translation
import periodic_atom_translation
from nonorthogonal_to_orthogonal_grid_interpolation import interpolate_using_linear_transformation, interp_linear_transform_orthogonal_z
from plot_dft_grid_data import plot_pot
from plot_atoms import plot_atoms

def find_atom(atom_types, atom_positions, atom_type_constraint, atom_pos_constraint=None):
    # Find all atoms matching the type
    type_matching_indices = [atom_index for atom_index, atom_type in enumerate(atom_types) if atom_type == atom_type_constraint]
    if len(type_matching_indices) == 1:
        return type_matching_indices[0]
    elif len(type_matching_indices) > 1:
        if atom_pos_constraint == 'max_z':
            return np.argmax(atom_positions[type_matching_indices, 2])
        elif atom_pos_constraint == 'min_z':
            return np.argmin(atom_positions[type_matching_indices, 2])
        else:
            print 'Found more than one atom matching the type but no supported position' + \
                'constraint was given. Returning the index of the first matching atom.'
            return type_matching_indices[0]
    else:
        sys.exit('Error: found no atoms with the matching type')


file_path = 'LOCPOT_2'
plot_filename = 'locpot_potential_2.pdf'
tip_atom_type = 8
grid_spacing = 0.08
plot_pos = [0, 0, 4.0]
plot_box_size = 20.0
atom_radius_scaling = 0.3
pot_xz_limits = [-8.0, 8.0]
pot_yz_limits = [-8.0, 8.0]
n_contour = 100

is_grid_orthogonal, is_z_orthogonal = locpot_import.check_grid_orthogonality(file_path)
if is_grid_orthogonal:
    print 'Reading LOCPOT file with orthogonal grid'
    atom_types, atom_positions, xs, ys, zs, pot_on_grid = locpot_import.read_locpot(file_path, orthogonal=True)
    print 'nx = {}, ny = {}, nz = {}'.format(xs.shape[0], ys.shape[0], zs.shape[0])
    anchor_atom = find_atom(atom_types, atom_positions, tip_atom_type, atom_pos_constraint='max_z')
    anchor_atom_pos = atom_positions[anchor_atom, :]
    print 'Translating the grid'
    xs, ys, zs, pot_on_grid = periodic_grid_translation.translate_grid(xs, ys, zs, pot_on_grid, anchor_atom_pos)
else:
    print 'Reading LOCPOT file non-orthogonal grid'
    atom_types, atom_positions, n_voxels, voxel_vectors, nonorthogonal_pot_grid = \
        locpot_import.read_locpot(file_path, orthogonal=False, return_unstructured_data=False)
    print 'nx_voxels = {}, ny_voxels = {}, nz_voxels = {}'.format(n_voxels[0], n_voxels[1], n_voxels[2])
    anchor_atom = find_atom(atom_types, atom_positions, tip_atom_type, atom_pos_constraint='max_z')
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

plt.figure(figsize=(8, 12))

ax1 = plt.subplot(2, 1, 1)
plot_pot(xs, ys, zs, pot_on_grid, plot_pos[1], slice_normal='y', limits=pot_xz_limits)
plot_atoms(atom_types, atom_positions, normal_direction='y', atom_radius_scaling=atom_radius_scaling)
plt.xlim(-plot_box_size/2, plot_box_size/2)
plt.ylim(-plot_box_size/2, plot_box_size/2)
ax1.set_aspect('equal')
plt.title(u'Electric potential (xz cut through tip)', size=16)
plt.xlabel(u'x (Å)', size=14)
plt.ylabel(u'z (Å)', size=14)

ax2 = plt.subplot(2, 1, 2)
plot_pot(xs, ys, zs, pot_on_grid, plot_pos[0], slice_normal='x', limits=pot_yz_limits)
plot_atoms(atom_types, atom_positions, normal_direction='x', atom_radius_scaling=atom_radius_scaling)
plt.xlim(-plot_box_size/2, plot_box_size/2)
plt.ylim(-plot_box_size/2, plot_box_size/2)
ax2.set_aspect('equal')
plt.title(u'Electric potential (yz cut through tip)', size=16)
plt.xlabel(u'y (Å)', size=14)
plt.ylabel(u'z (Å)', size=14)

plt.tight_layout(rect=(0, 0, 1, 0.95))
plt.savefig(plot_filename)
