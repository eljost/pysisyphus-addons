import numpy as np
import pytest
import time

from pysisyphus_addons.grid import grid
from pysisyphus.io.cube import get_grid
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.density_numba import eval_density as ref_eval_density


@pytest.mark.parametrize("num", [5, 10, 20, 40, 60])
def test_eval_density(num):
    """Compare Fortran implementation against numba implementation in pysisyphus."""
    wf_fn = "lib:orca_h2o_sto3g.json"
    wf = Wavefunction.from_file(wf_fn)
    grid3d, spacing, _ = get_grid(wf, num=num, margin=4.0)
    volume_element = np.prod(spacing)

    density = wf.P_tot
    reorder_c2s = wf.shells.reorder_c2s_coeffs
    # From external to pysisyphus
    density_cart = reorder_c2s.T @ density @ reorder_c2s

    fdur = time.time()
    dens_grid = grid.eval_density(wf.shells, grid3d, density_cart)
    fdur = time.time() - fdur
    N = dens_grid.sum() * volume_element
    print(f"{num=}, {dens_grid.shape=}, {N=: >12.6f}")

    # Compare against numba implementation
    shellstructs = wf.shells.as_numba_shellstructs()
    precontr = wf.shells.cart2sph_coeffs.T @ wf.shells.P_sph.T
    rho = np.zeros(len(grid3d))
    ndur = time.time()
    ref_eval_density(shellstructs, grid3d, density, precontr, rho)
    ndur = time.time() - ndur

    print(f"F took {fdur=: >12.6f}, {ndur=: >12.6f}")
    np.testing.assert_allclose(dens_grid, rho)


def test_densities():
    """Test multiple densities."""

    wf_fn = "lib:iao_ref_uhf.json"
    wf = Wavefunction.from_file(wf_fn)
    num = 80
    grid3d, spacing, _ = get_grid(wf, num=num, margin=4.0)
    volume_element = np.prod(spacing)

    Pa, Pb = wf.P
    P_tot = Pa + Pb
    P_spin = Pa - Pb
    reorder_c2s = wf.shells.reorder_c2s_coeffs
    # From external to pysisyphus
    density_tot = reorder_c2s.T @ P_tot @ reorder_c2s
    density_spin = reorder_c2s.T @ P_spin @ reorder_c2s
    densities = np.stack((density_tot, density_spin), axis=0)

    densities_grid = grid.eval_densities(wf.shells, grid3d, densities)
    N_refs = (9.972203, 0.0)
    for dens_grid, N_ref in zip(densities_grid, N_refs):
        N = dens_grid.sum() * volume_element
        print(f"{num=}, {dens_grid.shape=}, {N=: >12.6f}")
        assert N == pytest.approx(N_ref)
