# -*- coding: utf-8 -*-
#!/usr/bin/python

import pyximport; pyximport.install()
import sys
import argparse
import numpy as np
import cube_import
from nonorthogonal_to_orthogonal_grid_interpolation import interpolate_using_linear_transformation
from plot_dft_grid_data import plot_pot, plot_pot_yline, plot_pot_zline, plot_Ey_yline, plot_Ez_zline
import matplotlib.pyplot as plt

efield_relative_pos = 0.5
grid_spacing = 0.1
colormap = "sequential"
#plt.style.use("research_notes")

parser = argparse.ArgumentParser()
parser.add_argument("cube_file", help="The cube file that contains the electrostatic potential to be plotted.")
parser.add_argument("--min", type=float, default=0.0,
                    help="Minimum value of the potential after which rest are clipped.")
parser.add_argument("--max", type=float, default=0.0,
                    help="Maximum value of the potential after which rest are clipped.")
parser.add_argument("--nonorthorhombic", action="store_true",
                    help="Determines whether the unit cell represented by the cube file is nonorthorhombic.")
args = parser.parse_args()

if args.nonorthorhombic:
    print 'Reading cube file with non-orthogonal grid'
    atom_types, atom_pos, n_voxels, voxel_vectors, pot_on_nonortho_grid = cube_import.read_cube_nonorthogonal(args.cube_file, return_unstructured_data = False)
    print 'nx_voxels = {}, ny_voxels = {}, nz_voxels = {}'.format(n_voxels[0], n_voxels[1], n_voxels[2])
    print 'Interpolating to regular, orthogonal grid'
    new_origin = np.array([0.0, 0.0, 0.0])
    xs, ys, zs, pot_on_grid = interpolate_using_linear_transformation(pot_on_nonortho_grid,
                                            n_voxels, voxel_vectors, new_origin, grid_spacing)
    print 'nx = {}, ny = {}, nz = {}'.format(xs.shape[0], ys.shape[0], zs.shape[0])
else:
    print 'Reading cube file with orthogonal grid'
    atom_types, atom_pos, xs, ys, zs, pot_on_grid = cube_import.read_cube_orthogonal(args.cube_file)
    print 'nx = {}, ny = {}, nz = {}'.format(len(xs), len(ys), len(zs))

# Definition of Hartree potential
pot_on_grid = -pot_on_grid

efield_x = efield_relative_pos * (xs[1]-xs[0]) * xs.shape[0] + xs[0]
efield_y = efield_relative_pos * (ys[1]-ys[0]) * ys.shape[0] + ys[0]
efield_z = efield_relative_pos * (zs[1]-zs[0]) * zs.shape[0] + zs[0]

plt.figure(figsize=(11, 10))
if args.min == 0 and args.max == 0:
    plot_pot(xs, ys, zs, pot_on_grid, efield_z, slice_normal='z', colormap=colormap)
else:
    plot_pot(xs, ys, zs, pot_on_grid, efield_z, slice_normal='z', limits=[args.min, args.max], colormap=colormap)
plt.axes().set_aspect('equal')
plt.title(u'Electrostatic potential')
plt.xlabel(u'x (Å)')
plt.ylabel(u'y (Å)')

plt.figure(figsize=(10, 15))
plt.subplot(211)
ys_line, pot_line = plot_pot_yline(xs, ys, zs, pot_on_grid, np.array([efield_x, efield_z]), ys[0], ys[-1], n_points=450)
np.savetxt('pot_yline.txt', np.column_stack((ys_line, pot_line)))
plt.xlabel(u'y (Å)')
plt.ylabel(u'$\Phi$ (V)')
plt.subplot(212)
plot_Ey_yline(xs, ys, zs, pot_on_grid, np.array([efield_x, efield_z]), ys[0], ys[-1], n_points=500)
plt.xlabel(u'y (Å)')
plt.ylabel(u'$E_y$ (V/Å)')
plt.tight_layout()

plt.figure(figsize=(10, 15))
plt.subplot(211)
zs_line, pot_line = plot_pot_zline(xs, ys, zs, pot_on_grid, np.array([efield_x, efield_y]), zs[0], zs[-1], n_points=500)
np.savetxt('pot_zline.txt', np.column_stack((zs_line, pot_line)))
plt.xlabel(u'x (Å)')
plt.ylabel(u'$\Phi$ (V)')
plt.subplot(212)
plot_Ez_zline(xs, ys, zs, pot_on_grid, np.array([efield_x, efield_y]), zs[0], zs[-1], n_points=500)
plt.xlabel(u'x (Å)')
plt.ylabel(u'$E_x$ (V/Å)')
plt.tight_layout()

plt.show()
