import sys
import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy.distutils.system_info import default_include_dirs, default_lib_dirs

if sys.platform == 'win32':
    libraries = ["flint", "arb", "mpir", "mpfr", "pthreads"]
else:
    libraries = ["flint", "arb"]

default_include_dirs += [
    os.path.join(d, "flint") for d in default_include_dirs
]

ext_modules = [
    Extension(
        "flint", ["src/pyflint.pyx"],
        libraries=libraries,
        library_dirs=default_lib_dirs,
        include_dirs=default_include_dirs)
]

for e in ext_modules:
    e.cython_directives = {"embedsignature": True}

setup(
    name='python-flint',
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    description='bindings for FLINT',
    version='0.1.3',
    url='https://github.com/python-flint/python-flint',
    author='Fredrik Johansson',
    author_email='fredrik.johansson@gmail.com',
    license='BSD',
    classifiers=['Topic :: Scientific/Engineering :: Mathematics'])
