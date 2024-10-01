# Density evaluation on a grid, e.g., for cube generation

grid_mod = py.extension_module(
  '_grid',
  'grid.f90',
  'grid_wrapper.f90',
  '_grid.pyx',
  include_directories: np_inc,
  link_with: [libhermes],
  dependencies: [openmp_dep, py_dep],
  install: true,
  subdir: pkgdir / 'grid'
)

py.install_sources(
  '__init__.py',
  'grid.py',
  subdir: pkgdir / 'grid'
)
