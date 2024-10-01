cbd = meson.build_root()

test('Coulomb-Tensor for BeNa_2^2+',
  py,
  args: ['-m', 'pytest', '-s', '-v', files('test_coulomb_tensor.py')],
  env: ['PYTHONPATH=' + cbd + '/src'])
