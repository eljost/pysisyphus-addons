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
if buildtype == 'release'
  add_project_arguments([ '-Ofast', ], language: 'fortran')
endif

subdir('wavefunction')
