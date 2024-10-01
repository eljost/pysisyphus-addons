import numpy as np

from pysisyphus.wavefunction import Wavefunction
from pysisyphus_addons.diabatization import intor


def test_coulomb_tensor(this_dir):
    wf_fn = this_dir / "bena2_dcat.bson"
    wf = Wavefunction.from_file(wf_fn)
    shells = wf.shells
    densities_fn = this_dir / "bena2_dcat_densities_relaxed.npy"
    D = np.load(densities_fn)
    densities = np.stack((D[0, 0], D[0, 1], D[1, 1]), axis=0)
    coulomb_tensor = intor.contract_coulomb_densities_4d(shells, shells, densities)
    nstates = 2
    assert coulomb_tensor.shape == (nstates, nstates, nstates, nstates)
    coulomb_tensor_ref = np.array(
        [
            [
                [[344.547007476478, 0.00000000011], [0.00000000011, 344.802879544756]],
                [[0.00000000011, 0.014478048022], [0.014478048022, 0.000000000096]],
            ],
            [
                [[0.00000000011, 0.014478048022], [0.014478048022, 0.000000000096]],
                [
                    [344.802879544756, 0.000000000096],
                    [0.000000000096, 345.075649169291],
                ],
            ],
        ],
    )
    np.testing.assert_allclose(coulomb_tensor, coulomb_tensor_ref, atol=1e-12)
