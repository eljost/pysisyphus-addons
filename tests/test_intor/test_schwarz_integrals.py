import numpy as np
from pyscf import gto
import pytest

from pysisyphus_addons.intor import intor
from pysisyphus.wavefunction import Shells


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
    O    {La}
        1.0   1.0
        """
    )
    mol.basis = basis
    mol.build()

    shells = Shells.from_pyscf_mol(mol)

    # Quick comparison of overlap matrices
    S_ref = mol.intor("int1e_ovlp_sph")
    # The Overlap matrix is already correctly ordered
    S = shells.S_sph
    np.testing.assert_allclose(S, S_ref, atol=1e-14)

    int4c2e_ref = mol.intor("int2e_sph")
    # Reorder from PySCF to pysisyphus-order
    P = shells.P_sph
    int4c2e_ref_ro = np.einsum(
        "im,jn,ko,lp,ijkl->mnop", P, P, P, P, int4c2e_ref, optimize="greedy"
    )
    # Build matrix of Schwarz integrals
    naos = int4c2e_ref_ro.shape[0]
    schwarz_ref = np.zeros((naos, naos))
    for i in range(naos):
        for j in range(i, naos):
            schwarz_ref[i, j] = schwarz_ref[j, i] = int4c2e_ref_ro[i, j, i, j]

    # Calculate own matrix and compare. Basis functions are in pysisyphus' order.
    integrals = intor.int_schwarz(shells)
    np.testing.assert_allclose(integrals, schwarz_ref, atol=1e-14)
