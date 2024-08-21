project('python-flint', 'cython', 'c')

py = import('python').find_installation(pure: false)
dep_py = py.dependency()

cc = meson.get_compiler('c')

flint_ver_lower = '3.0' # >=3.0
flint_ver_upper = '3.2' # <3.2

gmp_dep = dependency('gmp')
mpfr_dep = dependency('mpfr')
flint_dep = dependency('flint', version: ['>=' + flint_ver_lower])

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
  if flint_dep.version().version_compare('>=' + flint_ver_upper)
    error('''

      Invalid Flint version:
      Version needed is: @0@ <= flint < @1@
      Version found is: flint == @2@

      By default, python-flint will only build against Flint versions that have
      been tested. If you are sure you want to use this version of Flint, you
      can disable this check with -Dflint_version_check=false.

      If building from the source directory using meson directly, you can do
      this with:

        meson setup build -Dflint_version_check=false

      If you are installing with pip, you can do this with:

        pip install --config-settings=setup-args="-Dflint_version_check=false" python-flint

      Other build frontends have similar options for passing this to meson.
      '''.format(flint_ver_lower, flint_ver_upper, flint_dep.version())
    )
  endif
endif

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
