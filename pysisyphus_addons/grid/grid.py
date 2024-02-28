from pysisyphus_addons.grid import _grid


def eval_density(shells, grid3d, density, **kwargs):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    return _grid.eval_density(
        bas_centers,
        bas_spec,
        bas_data,
        grid3d,
        density,
        **kwargs,
    )
