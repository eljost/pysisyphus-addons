import functools

from pysisyphus.wavefunction import Shells, Wavefunction

from pysisyphus_addons.std import _intor


@functools.singledispatch
def int_eri2c(shells: Shells, shells_aux: Shells):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    bas_centers_aux, bas_spec_aux, bas_data_aux = shells_aux.as_arrays(fortran=True)

    return _intor.int_eri2c(
        # AO basis
        bas_centers,
        bas_spec,
        bas_data,
        # Auxiliary basis
        bas_centers_aux,
        bas_spec_aux,
        bas_data_aux,
    )


@int_eri2c.register
def _(wf: Wavefunction):
    shells = wf.shells
    aux_basis_fn = "def2-universal-jfit.json"
    shells_aux = shells.from_basis(aux_basis_fn)
    return int_eri2c(shells, shells_aux)


def int_df1c(shells, shells_aux):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    bas_centers_aux, bas_spec_aux, bas_data_aux = shells_aux.as_arrays(fortran=True)

    return _intor.int_df1c(
        # AO basis
        bas_centers,
        bas_spec,
        bas_data,
        # Auxiliary basis
        bas_centers_aux,
        bas_spec_aux,
        bas_data_aux,
    )


def int_df2c_mo(shells, shells_aux, mo_tensor):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    bas_centers_aux, bas_spec_aux, bas_data_aux = shells_aux.as_arrays(fortran=True)

    return _intor.int_df2c_mo(
        # AO basis
        bas_centers,
        bas_spec,
        bas_data,
        # Auxiliary basis
        bas_centers_aux,
        bas_spec_aux,
        bas_data_aux,
        mo_tensor,
    )
