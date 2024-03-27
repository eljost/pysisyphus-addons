testdrive_dep = dependency('test-drive', fallback: ['test-drive', 'testdrive_dep'])

tests = [
  'erfinv',
  'screening',
]
test_srcs = files(
  'main.f90',
)
foreach t: tests
  test_srcs +=  files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [testdrive_dep],
  link_with: [libhermes],
)

foreach t: tests
  test(t, tester, args: t)
endforeach
