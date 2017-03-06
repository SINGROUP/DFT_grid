# -*- coding: utf-8 -*-
#!/usr/bin/python

import pyximport; pyximport.install()
import sys
import numpy as np
import cube_import
from plot_dft_grid_data import plot_pot, plot_pot_yline, plot_pot_zline, plot_Ey_yline, plot_Ez_zline
import matplotlib.pyplot as plt

slice_relative_pos = 0.5

if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Usage: python plot_macro_pot_xy.py <macro_pot_file>")

atom_types, atom_pos, xs, ys, zs, pot_on_grid = cube_import.read_cube_orthogonal(filename)
print 'nx = {}, ny = {}, nz = {}'.format(len(xs), len(ys), len(zs))

slice_x = slice_relative_pos * (xs[1]-xs[0]) * xs.shape[0] + xs[0]
slice_y = slice_relative_pos * (ys[1]-ys[0]) * ys.shape[0] + ys[0]
slice_z = slice_relative_pos * (zs[1]-zs[0]) * zs.shape[0] + zs[0]

plt.figure()
plot_pot(xs, ys, zs, pot_on_grid, slice_z, slice_normal='z')
plt.axes().set_aspect('equal')
plt.title(u'Electrostatic potential', size=16)
plt.xlabel(u'x (Å)', size=14)
plt.ylabel(u'y (Å)', size=14)

plt.figure()
plt.subplot(211)
ys_line, pot_line = plot_pot_yline(xs, ys, zs, pot_on_grid, np.array([slice_x, slice_z]), ys[0], ys[-1], n_points=450)
np.savetxt('pot_yline.txt', np.column_stack((ys_line, pot_line)))
plt.xlabel(u'y (Å)', size=14)
plt.ylabel(u'$\Phi$ (V)', size=14)
plt.subplot(212)
plot_Ey_yline(xs, ys, zs, pot_on_grid, np.array([slice_x, slice_z]), ys[0], ys[-1], n_points=500)
plt.xlabel(u'y (Å)', size=14)
plt.ylabel(u'$E_y$ (V/Å)', size=14)
plt.tight_layout()

plt.figure()
plt.subplot(211)
zs_line, pot_line = plot_pot_zline(xs, ys, zs, pot_on_grid, np.array([slice_x, slice_y]), zs[0], zs[-1], n_points=500)
np.savetxt('pot_zline.txt', np.column_stack((zs_line, pot_line)))
plt.xlabel(u'x (Å)', size=14)
plt.ylabel(u'$\Phi$ (V)', size=14)
plt.subplot(212)
plot_Ez_zline(xs, ys, zs, pot_on_grid, np.array([slice_x, slice_y]), zs[0], zs[-1], n_points=500)
plt.xlabel(u'x (Å)', size=14)
plt.ylabel(u'$E_x$ (V/Å)', size=14)
plt.tight_layout()

plt.show()
