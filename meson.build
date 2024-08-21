project('python-flint', 'cython', 'c')

#
# For the source release, we should by default fail for newer versions of Flint
# that are untested with a nice error message:
#
#   ../meson.build:10:12:
#   ERROR: Dependency lookup for flint with method 'pkgconfig' failed:
#   Invalid version, need 'flint' ['<3.2'] found '3.2.0'.
#
# We need an option to disable this though so that we can test newer versions
# of Flint. Also good to have an escape hatch for users since we don't know
# that future versions of Flint will not work.
#
if get_option('flint_version_check')
  flint_version_required = ['>=3.0', '<3.2']
else
  # We do know that Flint < 3.0 will not work so no point allowing it.
  flint_version_required = ['>=3.0']
endif

py = import('python').find_installation(pure: false)
dep_py = py.dependency()

cc = meson.get_compiler('c')

gmp_dep = dependency('gmp')
mpfr_dep = dependency('mpfr')
flint_dep = dependency('flint', version: flint_version_required)

add_project_arguments(
    '-X', 'embedsignature=True',
    '-X', 'emit_code_comments=True',
    language : 'cython'
)

if get_option('coverage')
   add_project_arguments('-X', 'linetrace=True', language : 'cython')
   add_project_arguments('-DCYTHON_TRACE=1', language : 'c')
endif

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
