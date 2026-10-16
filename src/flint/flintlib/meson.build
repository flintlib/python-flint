# Cython declarations for libflint. Needed on the install path so a
# third-party .pyx can `cimport flint.types.acb` and call acb_* on `.val`.
py.install_sources(
  [
    '__init__.py',
    'README.md',
  ],
  pure: false,
  subdir: 'flint/flintlib',
)

subdir('functions')
subdir('types')
