import math

import cython
import numpy as py_np
cimport numpy as cnp

cnp.import_array()


cdef extern:
    void f_eval_densities(
            # AO basis
            int nshells,
            int ndata,
            const int *bas_centers,
            const int *bas_spec,
            const double *bas_data,
            # Grid
            int npoints,
            const double *grid3d,
            # Densities
            int ndens,
            int nbfs,
            const double *densities,
            const double *grid_densities,
            int blk_size,
            double thresh,
            bint accumulate,
    )


def eval_densities(
    cython.int[:, ::1] bas_centers,
    cython.int[:, ::1] bas_spec,
    cython.double[::1] bas_data,
    cython.double[:, ::1] grid3d,
    cython.double[:, :, ::1] densities,
    cython.int blk_size = 100,
    cython.double thresh = 1e-8,
    cython.bint accumulate = False,
):
    ndens = densities.shape[0]
    nbfs = densities.shape[1]
    npoints = grid3d.shape[0]

    cdef:
        # AO basis
        # 2D
        cnp.ndarray f_bas_centers = py_np.asfortranarray(bas_centers)
        cnp.ndarray f_bas_spec = py_np.asfortranarray(bas_spec)
        # 1D
        cnp.ndarray f_bas_data = py_np.asfortranarray(bas_data)
        # Grid; transposed to shape (3, npoints)
        cnp.ndarray f_grid3d = py_np.asfortranarray(grid3d.T)
        # Densities
        cnp.ndarray f_densities = py_np.asfortranarray(densities)
        cnp.ndarray grid_densities = py_np.zeros((ndens, npoints), dtype="double", order="F")

    f_eval_densities(
        # AO basis
        len(bas_centers),
        len(bas_data),
        <int *> f_bas_centers.data,
        <int *> f_bas_spec.data,
        <double *> f_bas_data.data,
        # Grid
        npoints,
        <double *> f_grid3d.data,
        # Densities
        ndens,
        nbfs,
        <double *> f_densities.data,
        <double *> grid_densities.data,
        blk_size,
        thresh,
        accumulate,
    )
    return grid_densities
