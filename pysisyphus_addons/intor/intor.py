from pysisyphus_addons.intor import _intor


def int_schwarz_aux(shells):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    return _intor.int_schwarz_aux(bas_centers, bas_spec, bas_data)


def int_schwarz_bra(shells):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    return _intor.int_schwarz_bra(bas_centers, bas_spec, bas_data)
