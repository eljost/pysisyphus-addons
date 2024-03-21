openmp_dep = dependency('openmp')
blas_dep = dependency(get_option('blas'))
lapack_dep = dependency(get_option('lapack'))
# Maximum angular momenta for which integrals will be generated
lmax = get_option('lmax')
lauxmax = get_option('lauxmax')


py = import('python').find_installation(pure: false)
py_dep = py.dependency()

py_install_dir = py.get_install_dir()
message('py.install_dir', py_install_dir)
message('datadir', get_option('datadir'))

np_incdir= run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()
np_inc = include_directories(np_incdir)

py.install_sources(
  '__init__.py',
  subdir: pkgdir
)

buildtype = get_option('buildtype')
if meson.get_compiler('fortran').get_id() == 'gcc'
  if buildtype == 'release'
    add_project_arguments([ '-Ofast', ], language: 'fortran')
  elif buildtype == 'debug'
    add_project_arguments('-fbounds-check', language: 'fortran')
  endif
endif


subdir('gdma')

libhermes_srcs = []
subdir('common')
libhermes = library('libhermes', libhermes_srcs, 
  dependencies: [openmp_dep, blas_dep, lapack_dep],
)

subdir('intor')
subdir('diabatization')
subdir('grid')
