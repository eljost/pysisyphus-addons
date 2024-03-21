import math

import cython
import numpy as py_np
cimport numpy as cnp

cnp.import_array()


cdef extern:
    void f_int_schwarz(
            # AO basis
            int nshells,
            int ndata,
            const int *bas_centers,
            const int *bas_spec,
            const double *bas_data,
            # Integrals
            int nbfs,
            const double *integrals,
    )


def nbfs_from_bas_spec(cython.int[:, ::1] bas_spec):
    nbfs = 0
    for bs in bas_spec:
        nbfs += 2 * bs[1] + 1
    return nbfs


def int_schwarz(
    cython.int[:, ::1] bas_centers,
    cython.int[:, ::1] bas_spec,
    cython.double[::1] bas_data,
):
    nbfs = nbfs_from_bas_spec(bas_spec)
    print("There are", nbfs, "basis functions")

    cdef:
        # AO basis
        # 2D
        cnp.ndarray f_bas_centers = py_np.asfortranarray(bas_centers)
        cnp.ndarray f_bas_spec = py_np.asfortranarray(bas_spec)
        # 1D
        cnp.ndarray f_bas_data = py_np.asfortranarray(bas_data)
        # Integrals
        cnp.ndarray f_integrals = py_np.zeros(
            (nbfs, nbfs),
            dtype="double",
            order="F",
        )

    f_int_schwarz(
        # AO basis
        len(bas_centers),
        len(bas_data),
        <int *> f_bas_centers.data,
        <int *> f_bas_spec.data,
        <double *> f_bas_data.data,
        # Integrals
        nbfs,
        <double *> f_integrals.data,
    )
    return f_integrals
