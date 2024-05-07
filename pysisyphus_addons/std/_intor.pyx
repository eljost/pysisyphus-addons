import math

import cython
import numpy as py_np
cimport numpy as cnp

cnp.import_array()


cdef extern:
    void f_intor_eri2c(
            # AO basis
            int nshells,
            int ndata,
            const int *bas_centers,
            const int *bas_spec,
            const double *bas_data,
            # Auxiliary basis
            int nshells_aux,
            int ndata_aux,
            const int *bas_centers_aux,
            const int *bas_spec_aux,
            const double *bas_data_aux,
            int nbfs2,
            const double *eri2c_tensor,
    )
    void f_intor_df1c(
            # AO basis
            int nshells,
            int ndata,
            const int *bas_centers,
            const int *bas_spec,
            const double *bas_data,
            # Auxiliary basis
            int nshells_aux,
            int ndata_aux,
            const int *bas_centers_aux,
            const int *bas_spec_aux,
            const double *bas_data_aux,
            int nbfs,
            int nbfs_aux,
            const double *df1c_tensor,
    )


def int_eri2c(
    cython.int[:, ::1] bas_centers,
    cython.int[:, ::1] bas_spec,
    cython.double[::1] bas_data,
    cython.int[:, ::1] bas_centers_aux,
    cython.int[:, ::1] bas_spec_aux,
    cython.double[::1] bas_data_aux,
):
    nbfs = 0
    for shell in bas_spec:
        nbfs += (2 * shell[1] + 1)

    cdef:
        # AO basis
        # 2D
        cnp.ndarray f_bas_centers = py_np.asfortranarray(bas_centers)
        cnp.ndarray f_bas_spec = py_np.asfortranarray(bas_spec)
        # 1D
        cnp.ndarray f_bas_data = py_np.asfortranarray(bas_data)
        # Auxiliary AO basis
        # 2D
        cnp.ndarray f_bas_centers_aux = py_np.asfortranarray(bas_centers_aux)
        cnp.ndarray f_bas_spec_aux = py_np.asfortranarray(bas_spec_aux)
        # 1D
        cnp.ndarray f_bas_data_aux = py_np.asfortranarray(bas_data_aux)
        # Densities
        cnp.ndarray f_eri2c_tensor = py_np.zeros(
            (nbfs, nbfs),
            dtype="double",
            order="F",
        )

    f_intor_eri2c(
        # AO basis
        len(bas_centers),
        len(bas_data),
        <int *> f_bas_centers.data,
        <int *> f_bas_spec.data,
        <double *> f_bas_data.data,
        # Auxiliary basis
        len(bas_centers_aux),
        len(bas_data_aux),
        <int *> f_bas_centers_aux.data,
        <int *> f_bas_spec_aux.data,
        <double *> f_bas_data_aux.data,
        # ERI tensor
        nbfs,
        <double *> f_eri2c_tensor.data,
    )
    return f_eri2c_tensor


def int_df1c(
    cython.int[:, ::1] bas_centers,
    cython.int[:, ::1] bas_spec,
    cython.double[::1] bas_data,
    cython.int[:, ::1] bas_centers_aux,
    cython.int[:, ::1] bas_spec_aux,
    cython.double[::1] bas_data_aux,
):
    nbfs = 0
    for shell in bas_spec:
        nbfs += (2 * shell[1] + 1)
    nbfs_aux = 0
    for shell in bas_spec_aux:
        nbfs_aux += (2 * shell[1] + 1)

    cdef:
        # AO basis
        # 2D
        cnp.ndarray f_bas_centers = py_np.asfortranarray(bas_centers)
        cnp.ndarray f_bas_spec = py_np.asfortranarray(bas_spec)
        # 1D
        cnp.ndarray f_bas_data = py_np.asfortranarray(bas_data)
        # Auxiliary AO basis
        # 2D
        cnp.ndarray f_bas_centers_aux = py_np.asfortranarray(bas_centers_aux)
        cnp.ndarray f_bas_spec_aux = py_np.asfortranarray(bas_spec_aux)
        # 1D
        cnp.ndarray f_bas_data_aux = py_np.asfortranarray(bas_data_aux)
        # Densities
        cnp.ndarray f_df1c_tensor = py_np.zeros(
            (nbfs, nbfs_aux),
            dtype="double",
            order="F",
        )

    f_intor_df1c(
        # AO basis
        len(bas_centers),
        len(bas_data),
        <int *> f_bas_centers.data,
        <int *> f_bas_spec.data,
        <double *> f_bas_data.data,
        # Auxiliary basis
        len(bas_centers_aux),
        len(bas_data_aux),
        <int *> f_bas_centers_aux.data,
        <int *> f_bas_spec_aux.data,
        <double *> f_bas_data_aux.data,
        # DF tensor
        nbfs,
        nbfs_aux,
        <double *> f_df1c_tensor.data,
    )
    return f_df1c_tensor
