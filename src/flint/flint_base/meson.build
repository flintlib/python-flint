pkgdir = 'flint/flint_base'

pyfiles = [
  '__init__.py',
]

exts = [
  'flint_base',
  'flint_context',
]

py.install_sources(
  pyfiles,
  pure: false,
  subdir: pkgdir,
)

foreach ext : exts
  py.extension_module(
    ext,
    ext + '.pyx',
    dependencies: pyflint_deps,
    install: true,
    subdir: pkgdir,
  )
endforeach
