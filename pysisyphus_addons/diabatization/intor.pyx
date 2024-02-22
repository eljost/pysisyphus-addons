import math

import cython
import numpy as py_np
cimport numpy as cnp

cnp.import_array()


cdef extern:
    void f_contract_coulomb_densities_ao(
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
            # Densities
            int ndens,
            int nbfs,
            const double *densities,
            # Coulomb tensor
            int nstates,
            const double *coulomb_tensor,
    )


def contract_coulomb_densities_ao(
    cython.int[:, ::1] bas_centers,
    cython.int[:, ::1] bas_spec,
    cython.double[::1] bas_data,
    cython.int[:, ::1] bas_centers_aux,
    cython.int[:, ::1] bas_spec_aux,
    cython.double[::1] bas_data_aux,
    cython.double[:, :, ::1] densities,
):
    ndens = densities.shape[0]
    nbfs = densities.shape[1]
    # "nstates * (nstates + 1) / 2 = ndens" solved for nstates
    nstates = int(math.isqrt(8 * ndens + 1) / 2.0 - 0.5)
    assert nstates * (nstates + 1) // 2 == ndens

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
        cnp.ndarray f_densities = py_np.asfortranarray(densities)
        cnp.ndarray f_coulomb_tensor = py_np.zeros(
            (nstates, nstates, nstates, nstates),
            dtype="double",
            order="F",
        )
    

    f_contract_coulomb_densities_ao(
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
        # Densities
        ndens,
        nbfs,
        <double *> f_densities.data,
        # Coulomb tensor
        nstates,
        <double *> f_coulomb_tensor.data,
    )
    return f_coulomb_tensor
