import math

import cython
import numpy as py_np
cimport numpy as cnp

cnp.import_array()


cdef extern:
    void f_eval_density(
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
            int nbfs,
            const double *density,
            const double *grid_dens,
            int blk_size,
            double thresh,
            bint accumulate,
    )


def eval_density(
    cython.int[:, ::1] bas_centers,
    cython.int[:, ::1] bas_spec,
    cython.double[::1] bas_data,
    cython.double[:, ::1] grid3d,
    cython.double[:, ::1] density,
    cython.int blk_size = 100,
    cython.double thresh = 1e-8,
    cython.bint accumulate = False,
):
    nbfs = density.shape[1]
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
        cnp.ndarray f_density = py_np.asfortranarray(density)
        cnp.ndarray grid_dens = py_np.zeros(npoints, dtype="double")

    f_eval_density(
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
        nbfs,
        <double *> f_density.data,
        <double *> grid_dens.data,
        blk_size,
        thresh,
        accumulate,
    )
    return grid_dens
