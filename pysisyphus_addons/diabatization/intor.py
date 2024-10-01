from pysisyphus_addons.diabatization import _intor


def contract_coulomb_densities_2d(shells, shells_aux, densities):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    bas_centers_aux, bas_spec_aux, bas_data_aux = shells_aux.as_arrays(fortran=True)
    # bas_centers_aux, bas_spec_aux, bas_data_aux = shells.as_arrays(fortran=True)
    return _intor.contract_coulomb_densities_2d(
        # AO basis
        bas_centers,
        bas_spec,
        bas_data,
        # Auxiliary basis
        bas_centers_aux,
        bas_spec_aux,
        bas_data_aux,
        # Densities
        densities,
    )


def contract_coulomb_densities_4d(shells, shells_aux, densities):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    bas_centers_aux, bas_spec_aux, bas_data_aux = shells_aux.as_arrays(fortran=True)
    return _intor.contract_coulomb_densities_4d(
        # AO basis
        bas_centers,
        bas_spec,
        bas_data,
        # Auxiliary basis
        bas_centers_aux,
        bas_spec_aux,
        bas_data_aux,
        # Densities
        densities,
    )
