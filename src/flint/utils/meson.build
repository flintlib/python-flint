thisdir = 'flint/utils'

pyfiles = [
  '__init__.py',
  'flint_exceptions.py',
]

exts = []

py.install_sources(
  pyfiles,
  pure: false,
  subdir: thisdir,
)

foreach ext : exts
  py.extension_module(
    ext,
    ext + '.pyx',
    dependencies: pyflint_deps,
    install: true,
    subdir: thisdir,
  )
endforeach
