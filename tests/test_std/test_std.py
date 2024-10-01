from pysisyphus.wavefunction import Wavefunction, Shells

from pysisyphus_addons.std import intor

import numpy as np
from pyscf import df
import scipy as sp


np.set_printoptions(suppress=True, precision=4, linewidth=260)


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


def test_df2c_mo():
    fn = "lib:01_ch4_tzvp.json"
    fn = "/home/johannes/tmp/742_lih/01_lih_sto3g.bson"
    # fn = "/home/johannes/tmp/742_lih/02_h2_sto3g.bson"
    wf = Wavefunction.from_file(fn)
    shells = wf.shells
    mo_tensor, _ = wf.C
    # Bring basis functions in mo_tensor into pysisyphus' order
    # mo_tensor = shells.P_sph.T @ mo_tensor
    # mo_tensor = shells.P_sph @ mo_tensor
    mo_tensor = np.einsum("ri,ar", shells.P_sph, mo_tensor)
    # mo_tensor = mo_tensor[:, :2].copy()
    _, nmos = mo_tensor.shape
    # mo_tensor = np.eye(nmos)
    shells_ris = shells.to_ris_shells()
    df2c_mo = intor.int_df2c_mo(shells, shells_ris, mo_tensor)
    # df2c_mo = intor.int_df2c_mo(shells, shells, mo_tensor)

    # Calcualte DF integrals using pyscf
    mol = shells.to_pyscf_mol()
    auxmol = shells_ris.to_pyscf_mol()
    int3c = df.incore.aux_e2(mol, auxmol, intor="int3c2e_sph")
    print("boom")
    # int3c = mol.intor("int3c2e_sph")

    # Reorder PySCF integrals into pysisyphus-order
    pyscf_shells = Shells.from_pyscf_mol(mol)
    P_sph = pyscf_shells.P_sph
    # From pyscf to pysisyphus
    int3c = np.einsum("rsP,ri,sj->ijP", int3c, P_sph, P_sph, optimize="greedy")

    int3c_mo = np.einsum(
        "rsP,ri,sa->iaP", int3c, mo_tensor, mo_tensor, optimize="greedy"
    )
    int2c = auxmol.intor("int2c2e_sph")
    metric = np.linalg.pinv(int2c)
    eri_mo_ref = np.einsum("ijP,PQ,klQ", int3c_mo, metric, int3c_mo, optimize="greedy")
    L = np.linalg.cholesky(int2c)
    L_inv, info = sp.linalg.lapack.dtrtri(L, lower=1)
    print("metric")
    print(L_inv)
    # print(f"{info=}")
    B = int3c_mo @ L_inv.T
    # B = int3c_mo
    # B = int3c
    eri_mo_comp = np.einsum("ijP,klP->ijkl", B, B)
    np.testing.assert_allclose(eri_mo_comp, eri_mo_ref, atol=1e-10)
    ref_tensor = np.zeros_like(df2c_mo)
    for i in range(nmos):
        for a in range(i + 1):
            pack = (i + 1) * i // 2 + a
            print(f"{pack: >3d}: ({i: >2d}, {a: >2d})")
            ref_tensor[pack] = B[i, a]

    print("ref")
    print(ref_tensor)
    print("mine")
    print(df2c_mo)
    np.testing.assert_allclose(df2c_mo, ref_tensor, atol=1e-10)
