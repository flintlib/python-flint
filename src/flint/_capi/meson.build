thisdir = 'flint/_capi'

py.install_sources(
  '__init__.py',
  'cython.pxd',
  pure: false,
  subdir: thisdir,
)

py.extension_module(
  'cython',
  'cython.pyx',
  dependencies: pyflint_deps,
  install: true,
  subdir: thisdir,
  limited_api: limited_api_version,
)
