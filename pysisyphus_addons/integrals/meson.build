# 2-center-2-electron Coulomb integrals
int_2c2e = custom_target('int_2c2e',
  output: ['int_2c2e.f90', ],
  command: ['sympleints-graph', 'int2c2e',
            '--lmax', lmax.to_string(),
            '--lauxmax', lauxmax.to_string(),
            '--out-dir', pkgdir / 'integrals'])

# 3-center-2-electron Coulomb integrals
int_3c2e = custom_target('int_3c2e',
  output: ['int_3c2e.f90', ],
  command: ['sympleints-graph', 'int3c2e',
            '--lmax', lmax.to_string(),
            '--lauxmax', lauxmax.to_string(),
            '--out-dir', pkgdir / 'integrals'])

# Boys-function
libhermes_srcs += files(
  'boys_data.f90',
  'boys.f90',
)

# Generate density-fitting integrals
libhermes_srcs += [int_2c2e, int_3c2e]

libhermes_srcs += files(
  'init.f90',
  'intor.f90',
)
