# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys, os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from plot_locpot_efield_z import plot_efield_z

fig_collection_filename = 'CO_ir_tip_efields.pdf'
locpot_filename = 'LOCPOT'
base_folder = 'ir_co_tips'
locpot_files = [['14cluster', ''],
                ['30cluster', ''],
                ['cluster', ''],
                ['cluster100_surf', ''],
                ['cluster2100_surf', ''],
                ['co_14cluster', ''],
                ['co_30cluster', ''],
                ['co_clus2100surf', ''],
                ['co_cluster', ''],
                ['co_cluster100_surf', ''],
                ['co_l100_surf', ''],
                ['co_l111_surf', ''],
                ['co_lcluster', ''],
                ['co_lsurface_pyramid', ''],
                ['co_surface_pyramid', ''],
                ['l100_surface', ''],
                ['large_111_surface', ''],
                ['large_cluster', ''],
                ['lsurface_pyramid', ''],
                ['surface_pyramid', '']
                ]
grid_spacing = 0.08
plot_pos = [0, 0, 4.0]
plot_box_size = 12.0
atom_radius_scaling = 0.3
efield_xz_limits = [-0.3, 0.3]
efield_xy_limits = [-0.1, 0.1]

with PdfPages(fig_collection_filename) as pdf:
    for locpot_file in locpot_files:
        foldername = locpot_file[0]
        description = locpot_file[1]
        file_path = os.path.join(base_folder, foldername, locpot_filename)
        fig_file_path = os.path.join(base_folder, '{}.pdf'.format(foldername))
        print '- Plotting {} ({})'.format(description, foldername)
        
        plt.figure(figsize=(8, 12))
        plt.suptitle('{} ({})'.format(description, foldername), fontsize=16, fontweight='bold')
        plot_efield_z(file_path, grid_spacing, plot_pos, plot_box_size, efield_xy_limits, efield_xz_limits, atom_radius_scaling)
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.savefig(fig_file_path)
        pdf.savefig()
        plt.close()
        
        print 'Done\n'
