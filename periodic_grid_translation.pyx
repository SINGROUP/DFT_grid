import numpy as np
cimport numpy as np
cimport cython

DTYPE_FLOAT = np.float64
ctypedef np.float64_t DTYPE_FLOAT_t
DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

# Translates values on grid considering periodic boundary conditions
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def translate_grid(np.ndarray[DTYPE_FLOAT_t, ndim=1] xs,
                    np.ndarray[DTYPE_FLOAT_t, ndim=1] ys,
                    np.ndarray[DTYPE_FLOAT_t, ndim=1] zs,
                    np.ndarray[DTYPE_FLOAT_t, ndim=3] grid_data,
                    np.ndarray[DTYPE_FLOAT_t, ndim=1] new_center_pos):
    cdef int ix, iy, iz, nx, ny, nz, n, translated_ix, translated_iy, translated_iz
    cdef double dx, dy, dz
    cdef np.ndarray[DTYPE_INT_t, ndim=1] grid_center_point, grid_new_center, grid_translation, \
                                            translated_x_inds, translated_y_inds, translated_z_inds
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] translation, translated_xs, translated_ys, translated_zs
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=3] translated_grid_data
    
    nx = xs.shape[0]
    ny = ys.shape[0]
    nz = zs.shape[0]
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    dz = zs[1] - zs[0]
    grid_center_point = np.array([int(nx/2), int(ny/2), int(nz/2)])
    grid_new_center = np.array([int(new_center_pos[0]/dx), int(new_center_pos[1]/dy), int(new_center_pos[2]/dz)])
    grid_translation = grid_new_center - grid_center_point
    translation = np.array([grid_translation[0]*dx, grid_translation[1]*dy, grid_translation[2]*dz])
    
    translated_x_inds = np.zeros(nx, dtype=DTYPE_INT)
    translated_y_inds = np.zeros(ny, dtype=DTYPE_INT)
    translated_z_inds = np.zeros(nz, dtype=DTYPE_INT)
    for ix in range(nx):
        translated_ix = ix + grid_translation[0]
        if translated_ix < 0 or translated_ix >= nx:
            n = np.floor(float(translated_ix)/nx)
            translated_ix = translated_ix - n*nx
        translated_x_inds[ix] = translated_ix
    for iy in range(ny):
        translated_iy = iy + grid_translation[1]
        if translated_iy < 0 or translated_iy >= ny:
            n = np.floor(float(translated_iy)/ny)
            translated_iy = translated_iy - n*ny
        translated_y_inds[iy] = translated_iy
    for iz in range(nz):
        translated_iz = iz + grid_translation[2]
        if translated_iz < 0 or translated_iz >= nz:
            n = np.floor(float(translated_iz)/nz)
            translated_iz = translated_iz - n*nz
        translated_z_inds[iz] = translated_iz
    
    translated_grid_data = np.zeros((nx, ny, nz), dtype=DTYPE_FLOAT)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                translated_grid_data[ix, iy, iz] = grid_data[translated_x_inds[ix], translated_y_inds[iy], translated_z_inds[iz]]
    
    translated_xs = xs + translation[0]
    translated_ys = ys + translation[1]
    translated_zs = zs + translation[2]
    
    return translated_xs, translated_ys, translated_zs, translated_grid_data
