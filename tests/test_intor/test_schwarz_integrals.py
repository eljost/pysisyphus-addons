import numpy as np
from pyscf import gto
import pytest

from pysisyphus_addons.intor import intor
from pysisyphus.wavefunction import Shells


def get_starts(mol):
    starts = np.zeros(mol.nbas, dtype=int)
    for i in range(mol.nbas - 1):
        isize = 2 * mol.bas_angular(i) + 1
        starts[i + 1] = starts[i] + isize
    return starts


def get_bra_norms(mol, ord=None):
    int4c = mol.intor("int2e_sph")

    starts = get_starts(mol)
    shell_iter = range(mol.nbas)

    bra_norms = np.zeros((mol.nbas, mol.nbas))
    # Calculate quantities to estimate bra/principal basis shell pair
    for a in shell_iter:
        La = mol.bas_angular(a)
        sizea = 2 * La + 1
        slicea = slice(starts[a], starts[a] + sizea)
        for b in range(a, mol.nbas):
            Lb = mol.bas_angular(b)
            sizeb = 2 * Lb + 1
            sliceb = slice(starts[b], starts[b] + sizeb)
            int4c_slice = int4c[slicea, sliceb, slicea, sliceb]
            bra_norms[a, b] = bra_norms[b, a] = np.linalg.norm(int4c_slice, ord=ord)
    bra_norms = np.sqrt(bra_norms)
    return bra_norms


def get_aux_norms(mol, ord=None):
    int2c = mol.intor("int2c2e_sph")

    starts = get_starts(mol)
    shell_iter = range(mol.nbas)

    # Calculate quantities to estimate ket/aux. basis shell
    aux_norms = np.zeros(mol.nbas)
    for c in shell_iter:
        Lc = mol.bas_angular(c)
        sizec = 2 * Lc + 1
        slicec = slice(starts[c], starts[c] + sizec)
        ints = int2c[slicec, slicec]
        aux_norms[c] = np.linalg.norm(ints, ord=ord)
    aux_norms = np.sqrt(aux_norms)
    return aux_norms


@pytest.mark.parametrize("La", ("S", "P", "D", "F"))
@pytest.mark.parametrize("Lb", ("S", "P", "D", "F"))
def test_schwarz_integrals(La, Lb):
    mol = gto.Mole()
    mol.atom = (
        "O 0.0 -0.11081188 0.0; "
        "H 0.78397589 0.44324751 0.0; "
        "H -0.78397589    0.44324751    0.00000000;"
    )
    basis = gto.basis.parse(
        f"""
    H    {La}
         1.0              1.0
    H    {Lb}
          1.0               1.0
          0.5               0.5
    O    {La}
        1.0   1.0
        0.2   0.3
        """
    )
    mol.basis = basis
    mol.build()
    nshells = mol.nbas

    shells = Shells.from_pyscf_mol(mol)

    # Quick comparison of overlap matrices
    S_ref = mol.intor("int1e_ovlp_sph")
    # The Overlap matrix is already correctly ordered
    S = shells.S_sph
    np.testing.assert_allclose(S, S_ref, atol=1e-14)

    ref2d_bra = get_bra_norms(mol)
    # Compare only against lower triangle
    ref_bra = ref2d_bra[np.tril_indices(nshells)].flatten()

    # Calculate own matrix and compare. Basis functions are in pysisyphus' order.
    schwarz_bra = intor.int_schwarz_bra(shells)
    np.testing.assert_allclose(schwarz_bra, ref_bra, atol=1e-14)

    # Auxiliary Schwarz integrals
    ref_aux = get_aux_norms(mol)

    # Calculate own matrix and compare. Basis functions are in pysisyphus' order.
    schwarz_aux = intor.int_schwarz_aux(shells)
    np.testing.assert_allclose(schwarz_aux, ref_aux, atol=1e-14)
