from pysisyphus_addons.diabatization import intor


def contract_coulomb_densities_ao(shells, shells_aux, densities):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    bas_centers_aux, bas_spec_aux, bas_data_aux = shells_aux.as_arrays(fortran=True)
    coulomb_tensor = intor.contract_coulomb_densities_ao(
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
    return coulomb_tensor
