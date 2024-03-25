project('python-flint', 'cython', 'c')

py = import('python').find_installation(pure: false)
dep_py = py.dependency()

cc = meson.get_compiler('c')

gmp_dep = dependency('gmp')
mpfr_dep = dependency('mpfr')
flint_dep = dependency('flint')

flint_dep = dependency('flint')

if flint_dep.version().version_compare('<3.1')
  # This is a workaround for the fact that the flint.pc file in flint < 3.1
  # is missing -lflint. This is fixed in flint 3.1.0.
  flint_dep = cc.find_library('flint')
endif

pyflint_deps = [dep_py, gmp_dep, mpfr_dep, flint_dep]

subdir('src/flint')
