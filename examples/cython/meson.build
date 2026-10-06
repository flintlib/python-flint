project('python-flint-cython-example', 'cython', 'c', meson_version: '>=1.3')

py = import('python').find_installation(pure: false)
flint_dep = dependency('flint')

py.extension_module(
  '_flint_cython_example',
  '_flint_cython_example.pyx',
  dependencies: flint_dep,
  install: true,
)
