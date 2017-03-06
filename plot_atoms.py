# -*- coding: utf-8 -*-
#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.legend_handler import HandlerPatch
from ase.data import atomic_numbers, chemical_symbols, covalent_radii, vdw_radii
from ase.data.colors import jmol_colors

default_vdw_radius = 1.6

class HandlerCircle(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = Circle(xy=center, radius=height/2)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


def plot_atoms_height_xy(atom_types, atom_pos, z_offset=0.0, visible_atoms_set=None):
    if visible_atoms_set is None:
        visible_atoms_set = np.unique(atom_types)

    for visible_atom_type in visible_atoms_set:
        visible_atom_pos = atom_pos[(atom_types == visible_atom_type), :]
        plt.scatter(visible_atom_pos[:, 0], visible_atom_pos[:, 1],
                    s=visible_atom_type*15, c=visible_atom_pos[:, 2]-z_offset,
                    cmap=plt.cm.afmhot, alpha=0.4)
    
    cbar = plt.colorbar(pad=0.05, fraction=0.05, shrink=0.9)
    cbar.ax.set_ylabel(u'Height from Ir surface (Å)', size=14)


def plot_atoms(atom_types, atom_positions, visible_atoms_set=None, normal_direction='z', atom_plot_limits=None, atom_radius_scaling=1.0, atom_alpha=1.0):
    atom_min_height = 0
    atom_max_height = 0
    if normal_direction == 'x':
        atom_min_height = atom_positions[:, 0].min()
        atom_max_height = atom_positions[:, 0].max()
    elif normal_direction == 'y':
        atom_min_height = atom_positions[:, 1].min()
        atom_max_height = atom_positions[:, 1].max()
    elif normal_direction == 'z':
        atom_min_height = atom_positions[:, 2].min()
        atom_max_height = atom_positions[:, 2].max()
    else:
        raise Exception('normal_direction must be either x, y or z ({})'.format(normal_direction))
    
    ax = plt.gca()
    for ia in range(atom_types.shape[0]):
        atom_type = atom_types[ia]
        atom_pos = atom_positions[ia, :]
        vdw_radius = vdw_radii[atom_type]
        if np.isnan(vdw_radius):
            vdw_radius = default_vdw_radius
        if normal_direction == 'x':
            atom_normalized_height = (atom_pos[0] - atom_min_height) / (atom_max_height - atom_min_height)
            if atom_plot_limits is not None:
                if atom_pos[0] < atom_plot_limits[0] or atom_pos[0] > atom_plot_limits[1]:
                    continue
            atom_circle = Circle((atom_pos[1], atom_pos[2]), atom_radius_scaling*vdw_radius,
                                facecolor=jmol_colors[atom_type], zorder=atom_normalized_height+2, alpha=atom_alpha)
        elif normal_direction == 'y':
            atom_normalized_height = (atom_pos[1] - atom_min_height) / (atom_max_height - atom_min_height)
            if atom_plot_limits is not None:
                if atom_pos[1] < atom_plot_limits[0] or atom_pos[1] > atom_plot_limits[1]:
                    continue
            atom_circle = Circle((atom_pos[0], atom_pos[2]), atom_radius_scaling*vdw_radius,
                                facecolor=jmol_colors[atom_type], zorder=atom_normalized_height+2, alpha=atom_alpha)
        elif normal_direction == 'z':
            atom_normalized_height = (atom_pos[2] - atom_min_height) / (atom_max_height - atom_min_height)
            if atom_plot_limits is not None:
                if atom_pos[2] < atom_plot_limits[0] or atom_pos[2] > atom_plot_limits[1]:
                    continue
            atom_circle = Circle((atom_pos[0], atom_pos[1]), atom_radius_scaling*vdw_radius,
                                facecolor=jmol_colors[atom_type], zorder=atom_normalized_height+2, alpha=atom_alpha)
        ax.add_patch(atom_circle)
