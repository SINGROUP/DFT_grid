================
DFT grid atomic geo post proc
================
Description
-----------

Contains tools for post processing and visualizing grid based data from DFT calculations (charge density, electric potential) using Python. Data on a non-orthogonal grid can be efficiently interpolated into a regular orthogonal grid and the grid points can be translated inside the periodic cell. Contains many functions for plotting slices of different kinds of data represented on the grid. Functions for plotting projections of the atomic geometry are also included. Uses a separate `PyDFTFileIO <>`_ package to read the DFT output files containing the data grids.

Requirements
------------

- Python 2.7
- Cython (`http://cython.org/ <http://cython.org/>`_)
- NumPy (version >= 1.9)
- SciPy (version >= 0.14)
- matplotlib
- `PyDFTFileIO <>`_

For ``plot_atoms``:
- Atomic Simulation Environment (ASE), `https://wiki.fysik.dtu.dk/ase/ <https://wiki.fysik.dtu.dk/ase/>`_

Installation
------------

Put this directory containing the Python and Cython (.pyx) modules to your ``PYTHONPATH`` environment variable. The Cython modules should be compiled automatically at runtime as long as the Python scripts using them have the line
``import pyximport; pyximport.install()``
and you have Cython installed. See `http://docs.cython.org/en/latest/src/tutorial/cython_tutorial.html <http://docs.cython.org/en/latest/src/tutorial/cython_tutorial.html>`_ for more information.

Usage
-----

nonorthogonal_to_orthogonal_grid_interpolation:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Can be used to interpolate data on a non-orthogonal grid into a regular orthogonal grid for plotting purposes. The interpolation functions are computationally heavy and thus this module is written with Cython. See ``plot_locpot_efield_z`` and ``scripts/nonorthogonal_to_orthogonal_cube`` for example use cases.

plot_atoms:
^^^^^^^^^^^
Can be used to plot projection of an atomic geometry into a plane spanned by two Cartersian coordinate axes. ``plot_atoms`` different kinds of atoms with different colors and sizes. Supports projection to any Cartesian plane. ``plot_atoms_height_xy`` represents the height of the atoms with respect to the plane (hard coded to xy) by mapping the height to a color.

plot_atoms_mechafm:
^^^^^^^^^^^^^^^^^^^
The ``plot_atoms_pbc`` function is a more complicated version of the ``plot_atoms_height_xy`` supporting repetition of the atomic configuration accross periodic boundary conditions. It can plot atoms to any area by repeating the unit cell. Works well for plotting atom overlay to MechAFM (`https://github.com/SINGROUP/MechAFM <https://github.com/SINGROUP/MechAFM>`_) images, thus the name of the module. Should be combined with ``plot_atoms`` module.

plot_dft_grid_data:
^^^^^^^^^^^^^^^^^^^
Can be used to plot different kinds of data on grids that can be obtained as an output from DFT codes. Currently supports plotting of charge density, electric potential and electric field (calculated by taking numerical gradient of the potential). One can plot either a slice of the 3D grid or values along a line parallel to one of the coordinate axes. See the ``plot_epot_cube`` and ``plot_locpot_efield_z`` modules and the scripts in ``examples`` and ``scripts`` folders for usage examples.

Author
------
Juha Ritala (2016)
`jritala@gmail.com <mailto:jritala@gmail.com>`_

