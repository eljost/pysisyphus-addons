import math

import cython
import numpy as py_np
cimport numpy as cnp

cnp.import_array()


cdef extern:
    void f_int_schwarz_aux(
            # AO basis
            int nshells,
            int ndata,
            const int *bas_centers,
            const int *bas_spec,
            const double *bas_data,
            # Integrals
            const double *integrals,
    )

    void f_int_schwarz_bra(
            # AO basis
            int nshells,
            int ndata,
            const int *bas_centers,
            const int *bas_spec,
            const double *bas_data,
            # Integrals
            const double *integrals,
    )


"""
def nbfs_from_bas_spec(cython.int[:, ::1] bas_spec):
    nbfs = 0
    for bs in bas_spec:
        nbfs += 2 * bs[1] + 1
    return nbfs
"""


def int_schwarz_aux(
    cython.int[:, ::1] bas_centers,
    cython.int[:, ::1] bas_spec,
    cython.double[::1] bas_data,
):
    nshells = len(bas_centers)

    cdef:
        # AO basis
        # 2D
        cnp.ndarray f_bas_centers = py_np.asfortranarray(bas_centers)
        cnp.ndarray f_bas_spec = py_np.asfortranarray(bas_spec)
        # 1D
        cnp.ndarray f_bas_data = py_np.asfortranarray(bas_data)
        # Integrals
        cnp.ndarray f_integrals = py_np.zeros(nshells, dtype="double")

    f_int_schwarz_aux(
        # AO basis
        nshells,
        len(bas_data),
        <int *> f_bas_centers.data,
        <int *> f_bas_spec.data,
        <double *> f_bas_data.data,
        # Integrals
        <double *> f_integrals.data,
    )
    return f_integrals


def int_schwarz_bra(
    cython.int[:, ::1] bas_centers,
    cython.int[:, ::1] bas_spec,
    cython.double[::1] bas_data,
):
    nshells = len(bas_centers)
    nbra = nshells * (nshells + 1) // 2

    cdef:
        # AO basis
        # 2D
        cnp.ndarray f_bas_centers = py_np.asfortranarray(bas_centers)
        cnp.ndarray f_bas_spec = py_np.asfortranarray(bas_spec)
        # 1D
        cnp.ndarray f_bas_data = py_np.asfortranarray(bas_data)
        # Integrals
        cnp.ndarray f_integrals = py_np.zeros(nbra, dtype="double")

    f_int_schwarz_bra(
        # AO basis
        nshells,
        len(bas_data),
        <int *> f_bas_centers.data,
        <int *> f_bas_spec.data,
        <double *> f_bas_data.data,
        # Integrals
        <double *> f_integrals.data,
    )
    return f_integrals
