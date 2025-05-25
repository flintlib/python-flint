project(
  'python-flint',
  'cython',
  'c',
  meson_version : '>=1.1',
)
#
# The minimum versions are because we know that it will not work with earlier
# versions. The maximum versions are because python-flint was not tested
# against future versions that didn't exist at the time of release. In future
# if it seems like new releases do not always break the build of python-flint
# then we can consider not using a speculative upper version cap here.
#
flint_lower = '>=3.0'
flint_upper = '<3.3'
cython_lower = '>=3.0.11'
cython_upper = '<=3.1.0a1'

py = import('python').find_installation(pure: false)
dep_py = py.dependency()

cc = meson.get_compiler('c')
cy = meson.get_compiler('cython')

gmp_dep = dependency('gmp')
mpfr_dep = dependency('mpfr')
flint_dep = dependency('flint')

#
# For the source release, we should by default fail for new untested versions
# with a clear error message about the version mismatch.
#
# We need an option to disable this though so that we can test newer versions
# of Flint. Also good to have an escape hatch for users since we don't know
# that future versions of Flint will not work.
#
ver_message = '''

Invalid @0@ version:
Version needed is: @0@ @2@, @3@
Version found is: @0@ == @1@

By default, python-flint will only build against @0@ versions that have
been tested. If you are sure you want to use this version of @0@, you can
disable this check with -Dflint_version_check=false.

If building from the source directory using meson directly, you can do this
with:

  meson setup build -Dflint_version_check=false

If you are installing with pip, you can do this with:

  pip install --config-settings=setup-args="-Dflint_version_check=false" python-flint

Other build frontends have similar options for passing this to meson.

'''
if get_option('flint_version_check')
  if not (flint_dep.version().version_compare(flint_lower) and
          flint_dep.version().version_compare(flint_upper))
    error(ver_message.format('FLINT', flint_dep.version(), flint_lower, flint_upper))
  endif
  if not (cy.version().version_compare(cython_lower) and
          cy.version().version_compare(cython_upper))
    error(ver_message.format('Cython', cy.version(), cython_lower, cython_upper))
  endif
endif

# flint.pc was missing -lflint until Flint 3.1.0
if flint_dep.version().version_compare('<3.1')
  flint_dep = cc.find_library('flint')
  have_acb_theta = false
else
  have_acb_theta = true
endif

pyflint_deps = [dep_py, gmp_dep, mpfr_dep, flint_dep]

add_project_arguments(
    '-X', 'embedsignature=True',
    '-X', 'emit_code_comments=True',
    language : 'cython'
)

# Enable free-threading if Cython is new enough. The check should be
# >= 3.1.0a1 but meson gets confused by the a1 alpha release suffix.
# so we go with >= 3.1 (which will be correct once 3.1 is released).
cy = meson.get_compiler('cython')
if cy.version().version_compare('>=3.1')
  message('Enabling freethreading')
  add_project_arguments('-Xfreethreading_compatible=true', language : 'cython')
else
  message('Disabling freethreading')
endif

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

subdir('src/flint')
