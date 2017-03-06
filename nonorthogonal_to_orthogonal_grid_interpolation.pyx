import numpy as np
cimport numpy as np
cimport cython

DTYPE_FLOAT = np.float64
ctypedef np.float64_t DTYPE_FLOAT_t
DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def interpolate_using_linear_transformation(np.ndarray[DTYPE_FLOAT_t, ndim=3] nonorthogonal_grid_data,
                                            np.ndarray[DTYPE_INT_t, ndim=1] n_voxels,
                                            np.ndarray[DTYPE_FLOAT_t, ndim=2] voxel_vectors,
                                            np.ndarray[DTYPE_FLOAT_t, ndim=1] new_center_point, double new_grid_spacing):
    cdef int ix, iy, iz, nx, ny, nz, ix_sample, iy_sample, iz_sample, n
    cdef double d00, d01, d10, d11, d0, d1
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=2] basis_matrix_inv = np.zeros((3, 3), dtype=DTYPE_FLOAT)
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=2] corner_points = np.zeros((8, 3), dtype=DTYPE_FLOAT)
    cdef np.ndarray[DTYPE_INT_t, ndim=1] grid_point = np.zeros(3, dtype=DTYPE_INT)
    cdef np.ndarray[DTYPE_INT_t, ndim=1] sample_grid_point, n_cell_shifts
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] grid_center_point, translation, min_point, max_point, \
                                            position, pos_in_grid_basis, normalized_dist, normalized_rev_dist
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] xs, ys, zs
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=3] interpolation_sample = np.zeros((2, 2, 2), dtype=DTYPE_FLOAT)
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=3] grid_data
    
    # Basis transformation matrix from standard euclidean space to the space spanned by the voxel (unit cell) vectors
    basis_matrix_inv_t = np.linalg.inv(voxel_vectors).T
    
    # Create suitable orthogonal grid
    grid_center_point = (n_voxels[0]*voxel_vectors[0, :] + n_voxels[1]*voxel_vectors[1, :] + n_voxels[2]*voxel_vectors[2, :]) / 2
    translation = new_center_point - grid_center_point
    n = 0
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                corner_points[n, :] = ix*n_voxels[0]*voxel_vectors[0, :] + iy*n_voxels[1]*voxel_vectors[1, :] + \
                                        iz*n_voxels[2]*voxel_vectors[2, :] + translation
                n = n + 1
    
    min_point = np.amin(corner_points, axis=0)
    max_point = np.amax(corner_points, axis=0)
    
    nx = int((max_point[0]-min_point[0])/new_grid_spacing)
    xs = np.linspace(min_point[0], max_point[0], nx, dtype=DTYPE_FLOAT)
    ny = int((max_point[1]-min_point[1])/new_grid_spacing)
    ys = np.linspace(min_point[1], max_point[1], ny, dtype=DTYPE_FLOAT)
    nz = int((max_point[2]-min_point[2])/new_grid_spacing)
    zs = np.linspace(min_point[2], max_point[2], nz, dtype=DTYPE_FLOAT)
    
    # Interpolate the data in non-orthogonal grid at the grid points of the orthogonal grid
    grid_data = np.zeros((nx, ny, nz), dtype=DTYPE_FLOAT)
    for ix in range(nx):
        if ix % 10 == 0:
            print 'ix = {}/{}'.format(ix, nx)
        for iy in range(ny):
            for iz in range(nz):
                position = np.array([xs[ix], ys[iy], zs[iz]], dtype=DTYPE_FLOAT)
                pos_in_grid_basis = np.dot(basis_matrix_inv_t, position)
                grid_point = np.floor(pos_in_grid_basis).astype(DTYPE_INT)
                normalized_dist = pos_in_grid_basis - grid_point
                normalized_rev_dist = np.array([1.0, 1.0, 1.0]) - normalized_dist
                
                # Collect data samples around the interpolation point
                for ix_sample in range(2):
                    for iy_sample in range(2):
                        for iz_sample in range(2):
                            sample_grid_point = grid_point + np.array([ix_sample, iy_sample, iz_sample])
                            # Periodic boundary conditions
                            if sample_grid_point[0] < 0 or sample_grid_point[0] >= n_voxels[0]:
                                n = np.floor(float(sample_grid_point[0])/n_voxels[0])
                                sample_grid_point[0] = sample_grid_point[0] - n*n_voxels[0]
                            if sample_grid_point[1] < 0 or sample_grid_point[1] >= n_voxels[1]:
                                n = np.floor(float(sample_grid_point[1])/n_voxels[1])
                                sample_grid_point[1] = sample_grid_point[1] - n*n_voxels[1]
                            if sample_grid_point[2] < 0 or sample_grid_point[2] >= n_voxels[2]:
                                n = np.floor(float(sample_grid_point[2])/n_voxels[2])
                                sample_grid_point[2] = sample_grid_point[2] - n*n_voxels[2]
                            interpolation_sample[ix_sample, iy_sample, iz_sample] = \
                                nonorthogonal_grid_data[sample_grid_point[0], sample_grid_point[1], sample_grid_point[2]]
                
                d00 = normalized_rev_dist[0] * interpolation_sample[0, 0, 0] + normalized_dist[0] * interpolation_sample[1, 0, 0]
                d01 = normalized_rev_dist[0] * interpolation_sample[0, 0, 1] + normalized_dist[0] * interpolation_sample[1, 0, 1]
                d10 = normalized_rev_dist[0] * interpolation_sample[0, 1, 0] + normalized_dist[0] * interpolation_sample[1, 1, 0]
                d11 = normalized_rev_dist[0] * interpolation_sample[0, 1, 1] + normalized_dist[0] * interpolation_sample[1, 1, 1]
                d0 = normalized_rev_dist[1] * d00 + normalized_dist[1] * d10
                d1 = normalized_rev_dist[1] * d01 + normalized_dist[1] * d11
                grid_data[ix, iy, iz] = normalized_rev_dist[2] * d0 + normalized_dist[2] * d1
    
    return xs, ys, zs, grid_data


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def interp_linear_transform_orthogonal_z(np.ndarray[DTYPE_FLOAT_t, ndim=3] nonorthogonal_grid_data,
                                            np.ndarray[DTYPE_INT_t, ndim=1] n_voxels,
                                            np.ndarray[DTYPE_FLOAT_t, ndim=2] voxel_vectors,
                                            np.ndarray[DTYPE_FLOAT_t, ndim=1] new_center_point, double new_grid_spacing):
    cdef int ix, iy, iz, nx, ny, nz, ix_sample, iy_sample, n, z_grid_shift, iz_shifted, sample_grid_point_ix, sample_grid_point_iy
    cdef double dz
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=2] basis_xy_inv = np.zeros((2, 2), dtype=DTYPE_FLOAT)
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=2] corner_points = np.zeros((8, 3), dtype=DTYPE_FLOAT)
    cdef np.ndarray[DTYPE_INT_t, ndim=1] grid_point = np.zeros(2, dtype=DTYPE_INT)
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] grid_center_point, translation, min_point, max_point, \
                                            position, pos_in_grid_basis, normalized_dist, normalized_rev_dist
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] xs, ys, zs
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=2] interpolation_sample = np.zeros((2, 2), dtype=DTYPE_FLOAT)
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=3] grid_data
    
    basis_xy_inv_t = np.linalg.inv(voxel_vectors[0:2, 0:2]).T
    
    # Create suitable orthogonal grid
    grid_center_point = (n_voxels[0]*voxel_vectors[0, :] + n_voxels[1]*voxel_vectors[1, :] + n_voxels[2]*voxel_vectors[2, :]) / 2
    translation = new_center_point - grid_center_point
    n = 0
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):
                corner_points[n, :] = ix*n_voxels[0]*voxel_vectors[0, :] + iy*n_voxels[1]*voxel_vectors[1, :] + \
                                        iz*n_voxels[2]*voxel_vectors[2, :] + translation
                n = n + 1
    
    min_point = np.amin(corner_points, axis=0)
    max_point = np.amax(corner_points, axis=0)
    
    nx = int((max_point[0]-min_point[0])/new_grid_spacing)
    xs = np.linspace(min_point[0], max_point[0], nx, dtype=DTYPE_FLOAT)
    ny = int((max_point[1]-min_point[1])/new_grid_spacing)
    ys = np.linspace(min_point[1], max_point[1], ny, dtype=DTYPE_FLOAT)
    
    nz = n_voxels[2]
    dz = voxel_vectors[2, 2]
    z_grid_shift = int(round(min_point[2] / dz))
    zs = np.linspace(z_grid_shift*dz, (z_grid_shift + nz)*dz, nz, endpoint=False, dtype=DTYPE_FLOAT)
    
    # Interpolate the data in non-orthogonal grid at the grid points of the orthogonal grid
    grid_data = np.zeros((nx, ny, nz), dtype=DTYPE_FLOAT)
    for iz in range(nz):
        if iz % 100 == 0 and iz > 0:
            print 'iz = {}/{}'.format(iz, nz)
        iz_shifted = iz + z_grid_shift
        if iz_shifted < 0 or iz_shifted >= nz:
            n = np.floor(float(iz_shifted)/nz)
            iz_shifted = iz_shifted - n*nz
        for ix in range(nx):
            for iy in range(ny):
                position = np.array([xs[ix], ys[iy]], dtype=DTYPE_FLOAT)
                pos_in_grid_basis = np.dot(basis_xy_inv_t, position)
                grid_point = np.floor(pos_in_grid_basis).astype(DTYPE_INT)
                normalized_dist = pos_in_grid_basis - grid_point
                normalized_rev_dist = np.array([1.0, 1.0]) - normalized_dist
                
                # Collect data samples around the interpolation point
                for ix_sample in range(2):
                    sample_grid_point_ix = grid_point[0] + ix_sample
                    # Periodic boundary conditions
                    if sample_grid_point_ix < 0 or sample_grid_point_ix >= n_voxels[0]:
                        n = np.floor(float(sample_grid_point_ix)/n_voxels[0])
                        sample_grid_point_ix = sample_grid_point_ix - n*n_voxels[0]
                    for iy_sample in range(2):
                        sample_grid_point_iy = grid_point[1] + iy_sample
                        # Periodic boundary conditions
                        if sample_grid_point_iy < 0 or sample_grid_point_iy >= n_voxels[0]:
                            n = np.floor(float(sample_grid_point_iy)/n_voxels[0])
                            sample_grid_point_iy = sample_grid_point_iy - n*n_voxels[0]
                        interpolation_sample[ix_sample, iy_sample] = \
                            nonorthogonal_grid_data[sample_grid_point_ix, sample_grid_point_iy, iz_shifted]
                
                grid_data[ix, iy, iz] = normalized_rev_dist[0] * normalized_rev_dist[1] * interpolation_sample[0, 0] + \
                                        normalized_dist[0] * normalized_rev_dist[1] * interpolation_sample[1, 0] + \
                                        normalized_rev_dist[0] * normalized_dist[1] * interpolation_sample[0, 1] + \
                                        normalized_dist[0] * normalized_dist[1] * interpolation_sample[1, 1]
                                        
    return xs, ys, zs, grid_data
