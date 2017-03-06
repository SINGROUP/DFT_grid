import numpy as np
cimport numpy as np

DTYPE_FLOAT = np.float64
ctypedef np.float64_t DTYPE_FLOAT_t
DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

def translate_atoms(np.ndarray[DTYPE_FLOAT_t, ndim=2] atom_positions,
                    np.ndarray[DTYPE_INT_t, ndim=1] n_voxels,
                    np.ndarray[DTYPE_FLOAT_t, ndim=2] voxel_vectors,
                    np.ndarray[DTYPE_FLOAT_t, ndim=1] translation):
    cdef int ia, natoms
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=2] basis_matrix_inv = np.zeros((3, 3), dtype=DTYPE_FLOAT)
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] translated_atom_pos, pos_in_grid_basis, grid_origin
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=2] translated_atom_positions
    
    basis_matrix_inv_t = np.linalg.inv(voxel_vectors).T
    grid_origin = -(n_voxels[0]*voxel_vectors[0, :] + n_voxels[1]*voxel_vectors[1, :] + n_voxels[2]*voxel_vectors[2, :]) / 2
    
    natoms = atom_positions.shape[0]
    translated_atom_positions = np.zeros((natoms, 3), dtype=DTYPE_FLOAT)
    for ia in range(natoms):
        translated_atom_pos = atom_positions[ia, :] + translation - grid_origin
        pos_in_grid_basis = np.dot(basis_matrix_inv_t, translated_atom_pos)
        # Periodic boundary conditions
        if pos_in_grid_basis[0] < 0 or pos_in_grid_basis[0] >= n_voxels[0]:
            n = np.floor(float(pos_in_grid_basis[0])/n_voxels[0])
            pos_in_grid_basis[0] = pos_in_grid_basis[0] - n*n_voxels[0]
        if pos_in_grid_basis[1] < 0 or pos_in_grid_basis[1] >= n_voxels[1]:
            n = np.floor(float(pos_in_grid_basis[1])/n_voxels[1])
            pos_in_grid_basis[1] = pos_in_grid_basis[1] - n*n_voxels[1]
        if pos_in_grid_basis[2] < 0 or pos_in_grid_basis[2] >= n_voxels[2]:
            n = np.floor(float(pos_in_grid_basis[2])/n_voxels[2])
            pos_in_grid_basis[2] = pos_in_grid_basis[2] - n*n_voxels[2]
        
        translated_atom_positions[ia] = np.dot(voxel_vectors.T, pos_in_grid_basis) + grid_origin
    
    return translated_atom_positions
