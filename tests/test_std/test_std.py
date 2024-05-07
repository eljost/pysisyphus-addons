from pysisyphus.wavefunction import Wavefunction, Shells

from pysisyphus_addons.std import intor

import numpy as np


np.set_printoptions(suppress=True, precision=3, linewidth=260)


def test_eri2c():
    fn = "lib:01_ch4_tzvp.json"
    wf = Wavefunction.from_file(fn)
    shells = wf.shells
    # Calculate (aa|bb) using pysisyphus-addons
    eri2c = intor.int_eri2c(shells, shells)

    # Calculate reference integrals w/ PySCF using density fitting
    mol = shells.to_pyscf_mol()
    int3c = mol.intor("int3c2e_sph")

    # Reorder PySCF integrals into pysisyphus-order
    pyscf_shells = Shells.from_pyscf_mol(mol)
    P_sph = pyscf_shells.P_sph
    # Either resort pysisyphus-integrals to PySCF order or vice versa
    # eri2c = P_sph @ eri2c @ P_sph.T
    int3c = np.einsum("rsP,ri,sj->ijP", int3c, P_sph, P_sph, optimize="greedy")

    # Prune int3c (ab|P) to keep only (aa|P); this is faster than calculating
    # All ERIs and pruning them afterwards.
    ind = 0
    inds = list()
    sizes = list()
    for i in range(mol.nbas):
        inds.append(ind)
        L = mol.bas_angular(i)
        size = 2 * L + 1
        sizes.append(size)
        ind += size
    nshells = len(inds)
    tot_size = np.array(sizes).sum()
    naux = int3c.shape[2]
    int3c_flat = np.zeros((tot_size, naux))
    row_ind = 0
    for a in range(nshells):
        a_start = inds[a]
        a_size = sizes[a]
        a_slice = slice(a_start, a_start + a_size)
        slc = int3c[a_slice, a_slice]
        for i in range(sizes[a]):
            int3c_flat[row_ind + i] = slc[i, i]
        row_ind += a_size

    int2c = mol.intor("int2c2e_sph")
    metric = np.linalg.pinv(int2c)
    # Reconstruct ERI-tensor
    eri2c_ref = np.einsum("rP,PQ,sQ", int3c_flat, metric, int3c_flat)
    np.testing.assert_allclose(eri2c, eri2c_ref, atol=1e-12)
    mb = eri2c.nbytes / 1e6
    print(f"ERIs require {mb:.2f} MB")
