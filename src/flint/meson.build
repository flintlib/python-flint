thisdir = 'flint'

pyfiles = [
  '__init__.py',
  '__init__.pxd',
  '__init__.cython-31.pxd',
  'typing.py',
  'py.typed',
  'pyflint.pyi',
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

py.install_sources(
  '../python_flint/fmpz.h',
  pure: false,
  subdir: 'flint/include/python_flint',
)

foreach ext : exts
  py.extension_module(
    ext,
    ext + '.pyx',
    dependencies: pyflint_deps,
    install: true,
    subdir: thisdir,
    limited_api: limited_api_version,
  )
endforeach

foreach pkg : pkgs
  subdir(pkg)
endforeach
