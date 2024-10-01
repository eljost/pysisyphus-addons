import cython

import numpy as py_np
cimport numpy as cnp


cnp.import_array()


cdef extern:
    void f_eval_prim_density(
            int nprims,
            int npoints,
            int nbfs,
            const int *ls_inds,
            const double *primdata,
            const double *coords3d,
            const double *p,
            double switch,
            double *rho
    )


def eval_prim_density(
    cython.int[:, ::1] Ls_inds,
    cython.double[:, ::1] prim_data,
    cython.double[:, ::1] coords3d,
    cython.double[:, ::1] P_tot_cart,
    double switch,
    cython.double[::1] rho_pseudo
):
    cdef:
        cnp.ndarray f_Ls_inds = py_np.asfortranarray(Ls_inds.T)
        cnp.ndarray f_prim_data = py_np.asfortranarray(prim_data.T)
        cnp.ndarray f_coords3d = py_np.asfortranarray(coords3d.T)
        # Original order; density matrix is symmetric and rho_pseudo is 1d
        cnp.ndarray f_P_tot_cart = py_np.asfortranarray(P_tot_cart)
        cnp.ndarray f_rho_pseudo = py_np.asfortranarray(rho_pseudo)

    f_eval_prim_density(
        len(Ls_inds),  # nprims
        len(coords3d),  # npoints
        len(P_tot_cart),  # nbfs
        <int *> f_Ls_inds.data,
        <double *> f_prim_data.data,
        <double *> f_coords3d.data,
        <double *> f_P_tot_cart.data,
        switch,
        <double *> f_rho_pseudo.data,
    )
