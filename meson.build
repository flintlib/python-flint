project('python-flint', 'cython', 'c')

py = import('python').find_installation(pure: false)
dep_py = py.dependency()

cc = meson.get_compiler('c')

gmp_dep = dependency('gmp')
mpfr_dep = dependency('mpfr')
flint_dep = dependency('flint')

# Add rpaths for a local build of flint found via pkgconfig
#   https://github.com/mesonbuild/meson/issues/13046
if get_option('add_flint_rpath')
  flint_lib_dir = flint_dep.get_pkgconfig_variable('libdir')
  add_project_link_arguments(
    '-Wl,-rpath=' + flint_lib_dir,
    language: 'c',
  )
endif

# flint.pc was missing -lflint until Flint 3.1.0
if flint_dep.version().version_compare('<3.1')
  flint_dep = cc.find_library('flint')
  have_acb_theta = false
else
  have_acb_theta = true
endif

pyflint_deps = [dep_py, gmp_dep, mpfr_dep, flint_dep]

subdir('src/flint')
