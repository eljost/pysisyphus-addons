integral_srcs = []

# Loop over all integrals required for density fitting and screening.
foreach key : ['int2c2e', 'int3c2e', 'schwarz']
  srcs = run_command(
    'sympleints-graph',
    '--lmax', lmax.to_string(),
    '--lauxmax', lauxmax.to_string(),
    '--sources',
    key,
  ).stdout().strip().split()

  foreach fn : srcs
    message(fn)
  endforeach

  integral_srcs += custom_target(key + 'with sympleints-graph',
    output: srcs,
    command: [
      'sympleints-graph', key,
      '--lmax', lmax.to_string(),
      '--lauxmax', lauxmax.to_string(),
      '--out-dir', pkgdir / 'common' / 'integrals'
    ]
  )
endforeach

# Boys function sources
libhermes_srcs += files(
  'boys_data.f90',
  'boys.f90',
)

# Integral sources
libhermes_srcs += integral_srcs

# Initializer and driver sources
libhermes_srcs += files(
  'init.f90',
  'intor.f90',
)
