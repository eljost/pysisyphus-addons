from pysisyphus_addons.intor import _benchmark


def benchmark_int3c2e(shells, shells_aux):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    bas_centers_aux, bas_spec_aux, bas_data_aux = shells_aux.as_arrays(fortran=True)
    return _benchmark.benchmark_int3c2e(
        # AO basis
        bas_centers,
        bas_spec,
        bas_data,
        # Auxiliary basis
        bas_centers_aux,
        bas_spec_aux,
        bas_data_aux,
    )
