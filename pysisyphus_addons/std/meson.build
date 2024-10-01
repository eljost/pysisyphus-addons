#libstd = library('libstd', 'intor.f90', 'intor_wrapper.f90',
  #dependencies: [openmp_dep, blas_dep, lapack_dep],
  #link_with: [libhermes],
#)

std_intor_mod = py.extension_module(
  '_intor',
  'intor.f90',
  'intor_wrapper.f90',
  '_intor.pyx',
  include_directories: np_inc,
  link_with: [libhermes],
  dependencies: [py_dep],
  install: true,
  subdir: pkgdir / 'std'
)

py.install_sources(
  '__init__.py',
  'intor.py',
  subdir: pkgdir / 'std'
)
