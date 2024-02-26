# Integral-direct Coulomb-tensor build for Edmiston-Ruedenberg localization

dia_intor_mod = py.extension_module(
  '_intor',
  'intor.f90',
  'intor_wrapper.f90',
  '_intor.pyx',
  include_directories: np_inc,
  link_with: [libhermes],
  dependencies: [py_dep],
  install: true,
  subdir: pkgdir / 'diabatization'
)

py.install_sources(
  '__init__.py',
  'intor.py',
  subdir: pkgdir / 'diabatization'
)
