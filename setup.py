import sys
import os
from subprocess import check_call

from Cython.Distutils import build_ext
from Cython.Build import cythonize


if sys.platform == 'win32' and sys.version_info < (3, 12):
    from distutils.core import setup
    from distutils.extension import Extension
    from numpy.distutils.system_info import default_include_dirs, default_lib_dirs
    from distutils.sysconfig import get_config_vars
else:
    from setuptools import setup
    from setuptools.extension import Extension
    from sysconfig import get_config_vars
    default_include_dirs = []
    default_lib_dirs = []


libraries = ["flint"]


if sys.platform == 'win32':
    #
    # This is used in CI to build wheels with mingw64
    #
    if os.getenv('PYTHON_FLINT_MINGW64'):
        includedir = os.path.join(os.path.dirname(__file__), '.local', 'include')
        librarydir1 = os.path.join(os.path.dirname(__file__), '.local', 'bin')
        librarydir2 = os.path.join(os.path.dirname(__file__), '.local', 'lib')
        librarydirs = [librarydir1, librarydir2]
        default_include_dirs += [includedir]
        default_lib_dirs += librarydirs
        # Add gcc to the PATH in GitHub Actions when this setup.py is called by
        # cibuildwheel.
        os.environ['PATH'] += r';C:\msys64\mingw64\bin'
        libraries += ["mpfr", "gmp"]
    elif os.getenv('PYTHON_FLINT_MINGW64_TMP'):
        # This would be used to build under Windows against these libraries if
        # they have been installed somewhere other than .local
        libraries += ["mpfr", "gmp"]
    else:
        # For the MSVC toolchain link with mpir instead of gmp
        libraries += ["mpir", "mpfr", "pthreads"]
else:
    libraries = ["flint"]
    (opt,) = get_config_vars('OPT')
    os.environ['OPT'] = " ".join(flag for flag in opt.split() if flag != '-Wstrict-prototypes')


define_macros = []
compiler_directives = {
    'language_level': 3,
    'binding': False,
}


# Enable coverage tracing
if os.getenv('PYTHON_FLINT_COVERAGE'):
    define_macros.append(('CYTHON_TRACE', 1))
    compiler_directives['linetrace'] = True


packages = [
    'flint',
    'flint.flintlib',
    'flint.flint_base',
    'flint.types',
    'flint.functions',
    'flint.utils',
    'flint.test',
]


ext_files = [
    ("flint.pyflint", ["src/flint/pyflint.pyx"]),

    ("flint.flint_base.flint_base", ["src/flint/flint_base/flint_base.pyx"]),
    ("flint.flint_base.flint_context", ["src/flint/flint_base/flint_context.pyx"]),

    ("flint.types.fmpz", ["src/flint/types/fmpz.pyx"]),
    ("flint.types.fmpz_vec", ["src/flint/types/fmpz_vec.pyx"]),
    ("flint.types.fmpz_poly", ["src/flint/types/fmpz_poly.pyx"]),
    ("flint.types.fmpz_mpoly", ["src/flint/types/fmpz_mpoly.pyx"]),
    ("flint.types.fmpz_mat", ["src/flint/types/fmpz_mat.pyx"]),
    ("flint.types.fmpz_series", ["src/flint/types/fmpz_series.pyx"]),

    ("flint.types.fmpq", ["src/flint/types/fmpq.pyx"]),
    ("flint.types.fmpq_vec", ["src/flint/types/fmpq_vec.pyx"]),
    ("flint.types.fmpq_poly", ["src/flint/types/fmpq_poly.pyx"]),
    ("flint.types.fmpq_mat", ["src/flint/types/fmpq_mat.pyx"]),
    ("flint.types.fmpq_series", ["src/flint/types/fmpq_series.pyx"]),

    ("flint.types.nmod", ["src/flint/types/nmod.pyx"]),
    ("flint.types.nmod_poly", ["src/flint/types/nmod_poly.pyx"]),
    ("flint.types.nmod_mat", ["src/flint/types/nmod_mat.pyx"]),
    ("flint.types.nmod_series", ["src/flint/types/nmod_series.pyx"]),
    ("flint.types.nmod_mpoly", ["src/flint/types/nmod_mpoly.pyx"]),

    ("flint.types.fmpz_mod", ["src/flint/types/fmpz_mod.pyx"]),
    ("flint.types.fmpz_mod_poly", ["src/flint/types/fmpz_mod_poly.pyx"]),
    ("flint.types.fmpz_mod_mpoly", ["src/flint/types/fmpz_mod_mpoly.pyx"]),
    ("flint.types.fmpz_mod_mat", ["src/flint/types/fmpz_mod_mat.pyx"]),

    ("flint.types.fmpq_mpoly", ["src/flint/types/fmpq_mpoly.pyx"]),

    ("flint.types.fq_default", ["src/flint/types/fq_default.pyx"]),
    ("flint.types.fq_default_poly", ["src/flint/types/fq_default_poly.pyx"]),

    ("flint.types.arf", ["src/flint/types/arf.pyx"]),
    ("flint.types.arb", ["src/flint/types/arb.pyx"]),
    ("flint.types.arb_poly", ["src/flint/types/arb_poly.pyx"]),
    ("flint.types.arb_mat", ["src/flint/types/arb_mat.pyx"]),
    ("flint.types.arb_series", ["src/flint/types/arb_series.pyx"]),
    ("flint.types.acb", ["src/flint/types/acb.pyx"]),
    ("flint.types.acb_poly", ["src/flint/types/acb_poly.pyx"]),
    ("flint.types.acb_mat", ["src/flint/types/acb_mat.pyx"]),
    ("flint.types.acb_series", ["src/flint/types/acb_series.pyx"]),

    ("flint.types.dirichlet", ["src/flint/types/dirichlet.pyx"]),

    ("flint.functions.showgood", ["src/flint/functions/showgood.pyx"]),
]

ext_options = {
    "libraries" : libraries,
    "library_dirs" : default_lib_dirs,
    "include_dirs" : default_include_dirs,
    "define_macros" : define_macros,
}

ext_modules = []
for mod_name, src_files in ext_files:
    ext = Extension(mod_name, src_files, **ext_options)
    ext_modules.append(ext)

for e in ext_modules:
    e.cython_directives = {"embedsignature": True}


setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(ext_modules, compiler_directives=compiler_directives),
    packages=packages,
    package_dir={'': 'src'},
)
