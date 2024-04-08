import math

import cython
import numpy as py_np
cimport numpy as cnp

cnp.import_array()


cdef extern:
    void f_benchmark_int3c2e(
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
    )


def benchmark_int3c2e(
    cython.int[:, ::1] bas_centers,
    cython.int[:, ::1] bas_spec,
    cython.double[::1] bas_data,
    cython.int[:, ::1] bas_centers_aux,
    cython.int[:, ::1] bas_spec_aux,
    cython.double[::1] bas_data_aux,
):
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

    f_benchmark_int3c2e(
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
    )
