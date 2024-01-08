gdma_int_mod = py.extension_module(
  'prim_dens',
  'prim_dens.f90',
  'prim_dens.pyx',
  include_directories: np_inc,
  dependencies: [py_dep],
  install: true,
  subdir: pkgdir / 'wavefunction'
)

py.install_sources(
  '__init__.py',
  'gdma_int.py',
  subdir: pkgdir / 'wavefunction'
)
