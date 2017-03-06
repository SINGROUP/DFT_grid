# -*- coding: utf-8 -*-
#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

def plot_atoms_pbc(atom_types, atom_pos, unit_cell, area, xy_center, z_offset=None, visible_atoms_set=None, axes=None):
    n_atoms = atom_pos.shape[0]
    
    model_center = np.array([0.0, 0.0, 0.0])
    for position in atom_pos:
        model_center = model_center + position
    model_center = model_center/n_atoms
    
    if z_offset == None:
        z_offset = np.amin(atom_pos[:,2])
    
    model_origin = np.array([xy_center[0]-model_center[0], xy_center[1]-model_center[1], -z_offset])
    
    cell_a = unit_cell[:, 0]
    cell_b = unit_cell[:, 1]
    cell_a_len = np.sqrt(cell_a.dot(cell_a))
    cell_b_len = np.sqrt(cell_b.dot(cell_b))
    
    if visible_atoms_set is None:
        visible_atoms_set = np.unique(atom_types)
    
    if axes is None:
            axes = plt
    
    pbc_atom_pos = []
    pbc_atom_types = []
    #cell_a_shift_limits = [range(-int(np.amax(area)/cell_a_len)-1, int(np.amax(area)/cell_a_len)+2]
    #cell_b_shift_limits = [range(-int(np.amax(area)/cell_b_len)-1, int(np.amax(area)/cell_b_len)+2]
    cell_a_shift_limits = [0, 1]
    cell_b_shift_limits = [0, 1]
    for n_cell_a_shift in range(cell_a_shift_limits[0], cell_a_shift_limits[1]):
        for n_cell_b_shift in range(cell_b_shift_limits[0], cell_b_shift_limits[1]):
            for i_atom, position in enumerate(atom_pos):
                shifted_position = model_origin + position + n_cell_a_shift*cell_a + n_cell_b_shift*cell_b
                if (shifted_position[0] >= 0.0 and shifted_position[0] <= area[0] and
                    shifted_position[1] >= 0.0 and shifted_position[1] <= area[1]):
                    pbc_atom_pos.append(shifted_position)
                    pbc_atom_types.append(atom_types[i_atom])
    pbc_atom_pos = np.array(pbc_atom_pos)
    pbc_atom_types = np.array(pbc_atom_types)
    
    visible_atom_pos = pbc_atom_pos[np.in1d(pbc_atom_types, visible_atoms_set), :]
    visible_atom_types = pbc_atom_types[np.in1d(pbc_atom_types, visible_atoms_set)]

    S = axes.scatter(visible_atom_pos[:, 0], visible_atom_pos[:, 1],
                 s=visible_atom_types*14, c=visible_atom_pos[:, 2],
                 cmap=plt.cm.afmhot, alpha=0.2)

    return S

    #plt.xlim(0, area[0])
    #plt.ylim(0, area[1])
    #plt.axes().set_aspect('equal')
    #cbar = plt.colorbar()
    #plt.show()
