from typing import Tuple

import numpy as np

from pysisyphus.numint import MolGrid
from pysisyphus.wavefunction import Shells, Wavefunction
from pysisyphus.wavefunction.gdma_int import get_prim_data

from pysisyphus_addons.wavefunction.prim_dens import eval_prim_density


def get_diffuse_density(
    wf: Wavefunction, mol_grid: MolGrid, switch: float
) -> np.ndarray:
    # Convert pysisyphus shells to simple arrays that can be passed to Fortran.
    Ls_inds, primdata = get_prim_data(wf.shells)

    rho_pseudo = np.empty_like(mol_grid.weights)
    coords3d = mol_grid.xyz
    # Convert (spherical) density matrix to Cartesian matrix
    P_tot = wf.P_tot
    reorder_c2s = wf.shells.reorder_c2s_coeffs
    P_tot_cart = reorder_c2s.T @ P_tot @ reorder_c2s

    eval_prim_density(
        Ls_inds.astype(np.int32),
        primdata,
        coords3d,
        P_tot_cart,
        switch,
        rho_pseudo,
    )
    return rho_pseudo
