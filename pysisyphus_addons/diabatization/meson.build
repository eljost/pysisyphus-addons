# Integral-direct Coulomb-tensor build for Edmiston-Ruedenberg localization

dia_intor_mod = py.extension_module(
  'intor',
  'diabatization.f90',
  'diabatization_wrapper.f90',
  'intor.pyx',
  include_directories: np_inc,
  link_with: [libhermes],
  dependencies: [py_dep],
  install: true,
  subdir: pkgdir / 'diabatization'
)

py.install_sources(
  '__init__.py',
  'diabatization.py',
  subdir: pkgdir / 'diabatization'
)
