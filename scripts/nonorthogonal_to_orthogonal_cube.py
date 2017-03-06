# -*- coding: utf-8 -*-
#!/usr/bin/python

import pyximport; pyximport.install()
import numpy as np
import sys
import cube_import
import cube_file_writer as cube_writer
from nonorthogonal_to_orthogonal_grid_interpolation import interp_linear_transform_orthogonal_z

if len(sys.argv) == 4:
    filename_in = sys.argv[1]
    filename_out = sys.argv[2]
    grid_spacing = float(sys.argv[3])
else:
    sys.exit("Usage: python {} <cube file in> <cube file out> <grid spacing (in Å)>".format(sys.argv[0]))

atom_types, atom_pos, n_voxels, voxel_vectors, nonortho_datagrid = cube_import.read_cube_nonorthogonal(filename_in)
grid_center_point = (n_voxels[0]*voxel_vectors[0, :] + n_voxels[1]*voxel_vectors[1, :] + n_voxels[2]*voxel_vectors[2, :]) / 2
xs, ys, zs, data_on_grid = interp_linear_transform_orthogonal_z(nonortho_datagrid, n_voxels, voxel_vectors, grid_center_point, grid_spacing)
print 'nx = {}, ny = {}, nz = {}'.format(len(xs), len(ys), len(zs))
print 'xmin = {}, xmax = {}'.format(xs[0], xs[-1])
print 'test value = {}'.format(data_on_grid[len(xs)/2, len(ys)/2, len(zs)/2])
cube_writer.write_to_cube_with_atoms(filename_out, xs, ys, zs, data_on_grid, atom_types, atom_pos)
