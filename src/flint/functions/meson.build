pyfiles = [
  '__init__.py',
]

exts = [
  'showgood',
]

py.install_sources(
  pyfiles,
  pure: false,
  subdir: 'flint/functions',
)

foreach ext : exts
  py.extension_module(
    ext,
    ext + '.pyx',
    dependencies: pyflint_deps,
    install: true,
    subdir: 'flint/functions',
  )
endforeach
