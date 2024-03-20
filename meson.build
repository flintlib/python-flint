project('python-flint', 'cython', 'c')

py = import('python').find_installation(pure: false)
dep_py = py.dependency()

cc = meson.get_compiler('c')

gmp_dep = cc.find_library('gmp')
mpfr_dep = cc.find_library('mpfr')
flint_dep = cc.find_library('flint')

pyflint_deps = [dep_py, gmp_dep, mpfr_dep, flint_dep]

subdir('src/flint')
