# -*- coding: utf-8 -*-
#!/usr/bin/python

import pyximport; pyximport.install()
import sys
import numpy as np
import matplotlib.pyplot as plt
import cube_import
from plot_dft_grid_data import plot_Ey_yline

atom_ind_2 = None
if len(sys.argv) == 5:
    filename = sys.argv[1]
    atom_ind = int(sys.argv[2])-1
    y_min = float(sys.argv[3])
    y_max = float(sys.argv[4])
elif len(sys.argv) == 6:
    filename = sys.argv[1]
    atom_ind = int(sys.argv[2])-1
    atom_ind_2 = int(sys.argv[3])-1
    y_min = float(sys.argv[4])
    y_max = float(sys.argv[5])
else:
    sys.exit("Usage: python {} <cube_file> <atom_ind> [atom_ind_2] <y_min> <y_max>".format(sys.argv[0]))

atom_types, atom_pos, xs, ys, zs, pot_on_grid = cube_import.read_cube_orthogonal(filename)

point_x = atom_pos[atom_ind, 0]
point_z = atom_pos[atom_ind, 2]
y_min = y_min + atom_pos[atom_ind, 1]
y_max = y_max + atom_pos[atom_ind, 1]
plot_Ey_yline(xs, ys, zs, pot_on_grid, (point_x, point_z), y_min, y_max)

if atom_ind_2 is not None:
    point_x = atom_pos[atom_ind_2, 0]
    point_z = atom_pos[atom_ind_2, 2]
    plot_Ey_yline(xs, ys, zs, pot_on_grid, (point_x, point_z), y_min, y_max)
plt.show()
