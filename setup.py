from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("flint", ["pyflint.pyx"], libraries=["flint"])]

for e  in ext_modules:
    e.pyrex_directives  =  {"embedsignature":  True}

setup(
  name = 'flint',
  version = '0.1',
  cmdclass = {'build_ext':build_ext},
  ext_modules = ext_modules
)
