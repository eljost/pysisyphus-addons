from pysisyphus_addons.intor import _intor


def int_schwarz(shells):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    return _intor.int_schwarz(bas_centers, bas_spec, bas_data)
