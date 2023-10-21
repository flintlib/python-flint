Python-FLINT
============

Python extension module wrapping FLINT (Fast Library for Number Theory)
and Arb (arbitrary-precision ball arithmetic). Features:

* Integers, rationals, integers mod n
* Real and complex numbers with rigorous error tracking
* Polynomials, power series and matrices over all the above types
* Lots of mathematical functions

Documentation: http://fredrikj.net/python-flint/

Repository: https://github.com/flintlib/python-flint/

Author: Fredrik Johansson <fredrik.johansson@gmail.com>

Installation
------------

Currently python-flint supports CPython versions 3.9-3.12. For Windows (x86-64)
or OSX (x86-64 or arm64) or Linux (x86-64 `manylinux_2_17`) there are CPython
binary wheels for python-flint on PyPI. For these platforms python-flint can be
installed simply with `pip`

    pip install python-flint

Alternatively python-flint can be installed using `conda`

    conda install -c conda-forge python-flint

It is also possible to use python-flint with some PyPy versions. Binary wheels
are not provided for this on PyPI but can be installed with conda.

Build from source
-----------------

For other platforms or architectures installation needs to build from source.
First install FLINT 3.0.0. Note that as of python-flint 0.5.0 only this *exact*
version of FLINT will work. In future it is hoped that the version requirement
for python-flint will be FLINT >= 3.0.0 but at the time of writing 3.0.0 is the
newest version of FLINT that has only been released recently and is the only
version that is supported by python-flint.

See here for instructions on building FLINT:

* http://flintlib.org/

The latest release of Python-FLINT can then be built and installed using:

    pip install 'cython>=3' numpy wheel
    pip install --no-build-isolation python-flint

Python-FLINT can also be installed from a git checkout or a source archive
as follows:

    pip install 'cython>=3' numpy wheel
    pip install --no-build-isolation .

A script that builds and installs FLINT and python-flint that is tested on
Ubuntu can be found in the git repo here:

* https://github.com/flintlib/python-flint/blob/master/bin/pip_install_ubuntu.sh

See the documentation for further notes on building and installing
python-flint:

* https://fredrikj.net/python-flint/setup.html

Examples
-------------------------------------

Import Python-FLINT:

    >>> from flint import *

Number-theoretic functions:

    >>> fmpz(1000).partitions_p()
    24061467864032622473692149727991
    >>> fmpq.bernoulli(64)
    -106783830147866529886385444979142647942017/510

Polynomial arithmetic:

    >>> a = fmpz_poly([1,2,3]); b = fmpz_poly([2,3,4]); a.gcd(a * b)
    3*x^2 + 2*x + 1
    >>> a = fmpz_poly(list(range(10001))); b = fmpz_poly(list(range(10000))); a.gcd(a * b).degree()
    10000
    >>> x = fmpz_poly([0,1]); ((1-x**2)*(1+x**3)**3*(1+x+2*x)).factor()
    (-1, [(3*x + 1, 1), (x + (-1), 1), (x^2 + (-1)*x + 1, 3), (x + 1, 4)])

Matrix arithmetic:

    >>> fmpz_mat([[1,1],[1,0]]) ** 10
    [89, 55]
    [55, 34]
    >>> fmpq_mat.hilbert(10,10).det()
    1/46206893947914691316295628839036278726983680000000000

Numerical evaluation:

    >>> showgood(lambda: (arb.pi() * arb(163).sqrt()).exp() - 640320**3 - 744, dps=25)
    -7.499274028018143111206461e-13
    >>> showgood(lambda: (arb.pi() * 10**100 + arb(1)/1000).sin(), dps=25)
    0.0009999998333333416666664683

Numerical integration:

    >>> ctx.dps = 30
    >>> acb.integral(lambda x, _: (-x**2).exp(), -100, 100) ** 2
    [3.141592653589793238462643383 +/- 3.11e-28]

To do
-------------------------------------

* Write more tests and add missing docstrings
* Wrap missing flint types: finite fields, p-adic numbers, rational functions
* Vector or array types (maybe)
* Many convenience methods
* Write generic implementations of functions missing for specific FLINT types
* Proper handling of special values in various places (throwing Python
  exceptions instead of aborting, etc.)
* Various automatic conversions
* Conversions to and from external types (numpy, sage, sympy, mpmath, gmpy)
* Improved printing and string input/output
* IPython hooks (TeX pretty-printing etc.)

CHANGELOG
-------------

0.5.0

Important compatibility changes:

- gh-80, gh-94, gh-98: Switch from Flint 2.9 to Flint 3.
- gh-100: Supports Python 3.12 by using setuptools instead of numpy.distutils.

New features:

- gh-87: Adds `fmpz_mod_poly` type for polynomials over `fmpz_mod`.
- gh-85: Adds discrete logarithms to `fmpz_mod`.
- gh-83: Introduces the `fmpz_mod` type for multi-precision integer mods.

Bug fixes:

- gh-93: Fixes a bug with `pow(int, int, fmpz)` which previously gave incorrect
  results.
- gh-78, gh-79: minor fixes for the `nmod` type.

0.4.4

- gh-75, gh-77: finish bulk of the work in refactoring `python-flint` into
  submodules
- gh-72: The roots method of `arb_poly` is not supported. Use either the
  `complex_roots` method or `acb_roots(p).roots()` to get the old behaviour of
  returning the complex roots. The `roots` method on `fmpz_poly` and
  `fmpq_poly` now return integer and rational roots respectively. To access
  complex roots on these types, use the `complex_roots` method. For `acb_poly`,
  both `roots` and `complex_roots` behave the same
- gh-71: Include files in sdist and fix issue gh-70
- gh-67: Continue refactoring job to introduce submodules into `python-flint`

0.4.3

- gh-63: The `roots` method of `arb_poly`, and `nmod_poly` is no longer
  supported. Use `acb_roots(p).roots()` to get the old behaviour of returning
  the roots as `acb`. Note that the `roots` method of `fmpz_poly` and
  `fmpq_poly` currently returns the complex roots of the polynomial.
- gh-61: Start refactoring job to introduce submodules into `python-flint`

0.4.2

- gh-57: Adds manylinux wheels

0.4.1

- gh-47: Removes Linux wheels, updates instructions for building from source.

0.4.0

- gh-45: Adds wheels for Windows, OSX and manylinux but the Linux wheels are
  broken.

License
------------

Python-FLINT is licensed MIT. FLINT and Arb are LGPL v2.1+.
