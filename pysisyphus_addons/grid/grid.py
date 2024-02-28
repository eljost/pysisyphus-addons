from pysisyphus_addons.grid import _grid


def eval_density(shells, grid3d, density, **kwargs):
    """Wrapper that can be used when only one density is present."""
    assert density.ndim == 2
    densities = density[None, ...]
    densities_grid = eval_densities(shells, grid3d, densities, **kwargs)
    # Drop first dimension
    return densities_grid[0]


def eval_densities(shells, grid3d, densities, **kwargs):
    bas_centers, bas_spec, bas_data = shells.as_arrays(fortran=True)
    return _grid.eval_densities(
        bas_centers,
        bas_spec,
        bas_data,
        grid3d,
        densities,
        **kwargs,
    )
