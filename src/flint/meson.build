thisdir = 'flint'

pyfiles = [
  '__init__.py',
  'typing.py',
  'py.typed',
]

exts = [
  'pyflint',
]

pkgs = [
  'flint_base',
  'types',
  'functions',
  'utils',
  'test',
]

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

foreach pkg : pkgs
  subdir(pkg)
endforeach
