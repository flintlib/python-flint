from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("flint", ["pyflint.pyx"], libraries=["flint"])]

setup(
  name = 'flint',
  cmdclass = {'build_ext':build_ext},
  ext_modules = ext_modules
)
