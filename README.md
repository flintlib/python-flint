Python-FLINT
============

Python extension module wrapping FLINT (Fast Library for Number Theory)
and Arb (arbitrary-precision ball arithmetic). Features:

* Integers, rationals, integers mod n
* Real and complex numbers with rigorous error tracking
* Polynomials, power series and matrices over all the above types
* Lots of mathematical functions

Documentation: https://python-flint.readthedocs.io/en/latest/

Repository: https://github.com/flintlib/python-flint/

Author: Fredrik Johansson <fredrik.johansson@gmail.com>

Installation
------------

Currently python-flint supports CPython versions 3.11-3.13 and 3.13t
(free-threaded) and provides binaries on PyPI for the following platforms:

- Windows (x86-64)
- MacOS (x86-64, arm64)
- Linux (manylinux: x86-64, aarch64)

For these platforms python-flint can be installed simply with `pip`

    pip install python-flint

Alternatively python-flint can be installed using `conda`

    conda install -c conda-forge python-flint

Build from source
-----------------

For other platforms or architectures installation needs to build from source.
First install FLINT 3. Starting with python-flint 0.5.0 older versions of Flint
such as 2.9 are not supported any more. Note that as of Flint 3 Arb no longer
needs to be built separately as it is now merged into Flint.

As of e.g. Ubuntu 24.04 a new enough version of FLINT (at least version 3) can
be installed from the Ubuntu repos like

    sudo apt-get install libflint-dev

For older distros the version in the repos is too old and a newer version of
FLINT needs to be built. See here for instructions on building FLINT:

* http://flintlib.org/

A script that builds and installs FLINT on Ubuntu can be found here:

* https://github.com/flintlib/python-flint/blob/master/bin/install_flint_ubuntu.sh

The latest release of Python-FLINT can then be built from source and installed
using:

    pip install --no-binary python-flint python-flint

Python-FLINT can also be installed from a git checkout or a source archive
as follows:

    pip install .

See the documentation for further notes on building and installing
python-flint:

* https://python-flint.readthedocs.io/en/latest/build.html
* https://python-flint.readthedocs.io/en/latest/install.html

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
* Wrap missing flint types: matrices over finite fields, p-adic numbers, rational functions
* Build on the preliminary interface to FLINT's generic (gr) types.
* Make a nicer interface like `ZZ(1)` etc rather than `fmpz_poly([1, 2])`.
* Vector or array types (maybe)
* Many convenience methods
* Write generic implementations of functions missing for specific FLINT types
* Proper handling of special values in various places (throwing Python
  exceptions instead of aborting, etc.)
* Various automatic conversions
* Conversions to and from external types (numpy, sage, sympy, mpmath, gmpy)
* Improved printing and string input/output
* IPython hooks (TeX pretty-printing etc.)

Compatibility table
-------------------

Generally each release of python-flint will be compatible with a range of
Python versions as described in [SPEC
0](https://scientific-python.org/specs/spec-0000/). Since python-flint 0.5.0
the minimum supported Flint version is `3.0` and each release of python-flint
supports all versions of Flint `>=3.0` available at the time of release.

Compatible versions (note that 0.7.0 is not yet released):

| python-flint | Release date  | CPython     | FLINT      | Cython           |
|--------------|---------------|-------------|------------|------------------|
| `0.7.0`      | 16th Mar 2025 | `3.11-3.13` | `3.0-3.2`  | `3.0.11-3.1.0a1` |
| `0.6.0`      |  1st Feb 2024 | `3.9-3.12`  | `3.0` only | `3.0` only       |

As of python-flint 0.7.0, CPython 3.13 [PEP
703](https://peps.python.org/pep-0703/) free-threaded (no-GIL) builds of
python-flint are provided. In the the free-threaded build, mutating matrices or
polynomials from multiple threads can lead to memory corruption. Provided
matrices or polynomials are not mutated when shared across threads there are no
known issues with the free-threaded build but these should be considered
experimental.

CHANGELOG
=========

Next release (0.8.0)...
-----------------------

0.7.0
-----

Contributors (0.7.0):

- Jake Moss (JM)
- Giacomo Pope (GP)
- Joris Roos (JR)
- Edgar Costa (EC)
- Frédéric Chapoton (FC)
- Oscar Benjamin (OB)
- Tom Hubrecht (TH)

Highlights (0.7.0):

- [gh-270](https://github.com/flintlib/python-flint/pull/270),
  PyPI packages are now built with FLINT 3.2.0 (previously
  3.0.1 was used). All versions from FLINT 3.0.0 to FLINT 3.2.0
  are compatible with python-flint but some features require
  newer FLINT versions and the PyPI packages now use FLINT 3.2.0.
- [gh-97](https://github.com/flintlib/python-flint/pull/97),
  [gh-182](https://github.com/flintlib/python-flint/pull/180):
  Add `fq_default` and `fq_default_poly` for finite fields and
  univariate polynomials over finite fields. This exposes all
  of the different implementations of finite fields (`fq_zech`,
  `fq_nmod` etc) via the `fq_default` interface. (GP)
- [gh-132](https://github.com/flintlib/python-flint/pull/132),
  [gh-164](https://github.com/flintlib/python-flint/pull/164),
  [gh-190](https://github.com/flintlib/python-flint/pull/190),
  [gh-191](https://github.com/flintlib/python-flint/pull/191):
  [gh-192](https://github.com/flintlib/python-flint/pull/192):
  [gh-216](https://github.com/flintlib/python-flint/pull/216):
  [gh-225](https://github.com/flintlib/python-flint/pull/225):
  [gh-228](https://github.com/flintlib/python-flint/pull/228):
  Add `fmpz_mpoly`, `fmpq_mpoly`, `nmod_poly` and `fmpz_mod_poly`
  types for multivariate polynomials with integer, rational or
  integers mod `n` coefficients. (JM)
- [gh-142](https://github.com/flintlib/python-flint/pull/142)
  Add `acb_theta` module for the numerical evaluation of [theta
  functions](https://flintlib.org/doc/acb_theta.html) (only
  available for `Flint >= 3.1`). (EC)
- [gh-218](https://github.com/flintlib/python-flint/pull/218)
  [gh-254](https://github.com/flintlib/python-flint/pull/254)
  [gh-255](https://github.com/flintlib/python-flint/pull/255)
  An experimental interface for FLINT's generic rings has been
  added. This provides access to many of FLINT's types that
  are not yet wrapped by python-flint such as Gaussian integer,
  number fields, qqbar, calcium, as well as both univariate and
  multivariate polynomials and series over these rings (no
  matrices yet though). (OB and TH)
- [gh-129](https://github.com/flintlib/python-flint/pull/129)
  [gh-208](https://github.com/flintlib/python-flint/pull/208)
  Use meson/meson-python instead of setuptools as the build system
  for parallel builds and better detection of build and dependency
  requirements. (OB)
- [gh-201](https://github.com/flintlib/python-flint/pull/201)
  [gh-202](https://github.com/flintlib/python-flint/pull/202)
  The documentation has been updated and is now at
  [readthedocs](https://python-flint.readthedocs.io/en/latest/).
  (OB)
  [gh-235](https://github.com/flintlib/python-flint/pull/235)
  Nightly wheels for python-flint can now be installed from the
  [Anaconda Scientific Python Nightly Wheels index]
  (https://anaconda.org/scientific-python-nightly-wheels/python-flint).
  [gh-259](https://github.com/flintlib/python-flint/pull/259)
  Add PyPI wheels for Linux aarch64 (Linux on ARM CPU). (OB)

Compatibility break (0.7.0):

- [gh-189](https://github.com/flintlib/python-flint/pull/189)
  As of python-flint 0.7.0 `fmpq_poly.factor()` now returns
  primitive rather than monic factors i.e. `2*x + 1` rather than
  `x + 1/2`. This ensures consistency between all poly types
  including between `fmpq_poly` and `fmpq_mpoly`. (OB)

Other changes (0.7.0):

- [gh-269](https://github.com/flintlib/python-flint/pull/269)
  All univariate and multivariate polynomial types have
  `is_zero`, `is_one` and `is_constant` methods. All polynomial
  types now consistently handle negative powers where possible.
- [gh-261](https://github.com/flintlib/python-flint/pull/261)
  Add `fmpz_mat.fflu` for fraction-free LU decomposition of
  an integer matrix.
- [gh-251](https://github.com/flintlib/python-flint/pull/251)
  Add mpmath-style precision context managers for arb
  `extraprec`, `extradps`, `workprec` and `workdps`. (TH)
- [gh-250](https://github.com/flintlib/python-flint/pull/250)
  Add `fmpq.gcd()` method.
- [gh-215](https://github.com/flintlib/python-flint/pull/215)
  [gh-219](https://github.com/flintlib/python-flint/pull/219)
  The FLINT binding declarations are now fully generated
  automatically from the FLINT docs. (OB)
- [gh-203](https://github.com/flintlib/python-flint/pull/203)
  [gh-204](https://github.com/flintlib/python-flint/pull/204)
  [gh-205](https://github.com/flintlib/python-flint/pull/205)
  [gh-206](https://github.com/flintlib/python-flint/pull/206)
  [gh-207](https://github.com/flintlib/python-flint/pull/207)
  [gh-211](https://github.com/flintlib/python-flint/pull/211)
  [gh-212](https://github.com/flintlib/python-flint/pull/212)
  [gh-271](https://github.com/flintlib/python-flint/pull/271)
  Various linting fixes and codebase improvements (FC and GP).
- [gh-189](https://github.com/flintlib/python-flint/pull/189)
  All scalar and poly types now have `sqrt`. All poly types now
  have `factor_squarefree` and `leading_coefficient` methods.
  Exception types raised in a number of places were changed to
  `DomainError` for better consistency. (OB)
- [gh-196](https://github.com/flintlib/python-flint/pull/196)
  Supported Python versions are 3.10-3.13 (3.9 dropped). CI
  Testing added for 3.13 free-threaded CPython.
- [gh-194](https://github.com/flintlib/python-flint/pull/194)
  Add version checking for build requirements. (OB)
- [gh-180](https://github.com/flintlib/python-flint/pull/180)
  Add `equal_trunc`, `add_trunc`, `sub_trunc`, `mul_low`,
  `mul_mod` and `pow_trunc` methods to `fmpz_mod_poly`. (GP)
- [gh-177](https://github.com/flintlib/python-flint/pull/177)
  Remove old Py2 code for compatibility with Cython 3.1. (OB)
- [gh-176](https://github.com/flintlib/python-flint/pull/176)
  Fix the error messages from `fmpq` constructor. (OB)
- [gh-174](https://github.com/flintlib/python-flint/pull/174)
  Add `pow_mod` and `compose_mod` methods to `nmod_poly` and
  `fmpz_mod_poly`. Also add some missing methods to `nmod_poly`
  that other poly types already have. (GP)
- [gh-172](https://github.com/flintlib/python-flint/pull/172)
  Add `fmpz_is_square`. (JR)
- [gh-168](https://github.com/flintlib/python-flint/pull/168)
  Make comparisons consistent between different types. Add
  `is_one` and `is_zero` for all poly types. (OB)
- [gh-161](https://github.com/flintlib/python-flint/pull/161)
  Add `acb.lerch_phi` to compute the Lerch transcendent. (OB)
- [gh-160](https://github.com/flintlib/python-flint/pull/160)
  Add `bits` to `arb` and `acb`, add `log_base` to `arb`. (JR)
- [gh-148](https://github.com/flintlib/python-flint/pull/148)
  Remove debug symbols to make smaller Linux binaries. (OB)
- [gh-144](https://github.com/flintlib/python-flint/pull/144)
  Add `rel_one_accuracy_bits` to `arb` and `acb`. (EC)
- [gh-137](https://github.com/flintlib/python-flint/pull/137)
  Add `erfinv` and `erfcinv` for `arb`. (JR)
- [gh-119](https://github.com/flintlib/python-flint/pull/119)
  Add compatibility with Flint 3.1. (OB)

0.6.0
-----

- [gh-112](https://github.com/flintlib/python-flint/issues/112),
  [gh-111](https://github.com/flintlib/python-flint/issues/111),
  [gh-110](https://github.com/flintlib/python-flint/issues/110),
  [gh-108](https://github.com/flintlib/python-flint/issues/108):
  Add pyproject.toml and build dependencies. This means that
  python-flint can be built from source without
  `--no-build-isolation`.
- [gh-109](https://github.com/flintlib/python-flint/issues/109):
  Use exact division for non-field domains. Now `fmpz(6)/fmpz(3)`
  returns an exact result `fmpz(2)` or raises an error if an exact
  result is not possible. Similar changes for `fmpz_poly/fmpz`,
  `fmpz_mat/fmpz`, and for polynomial division with `fmpz_poly`,
  `fmpq_poly`, `nmod_poly` and `fmpz_mod_poly`.
- [gh-106](https://github.com/flintlib/python-flint/issues/106):
  Add `fmpz_mod_mat` for matrices of integers mod `n` where `n` is
  larger than word sized.
- [gh-104](https://github.com/flintlib/python-flint/issues/104):
  Bump Flint from 3.0.0 to 3.0.1

0.5.0
-----

Important compatibility changes:

- [gh-80](https://github.com/flintlib/python-flint/issues/80),
  [gh-94](https://github.com/flintlib/python-flint/issues/94),
  [gh-98](https://github.com/flintlib/python-flint/issues/98):
  Switch from Flint 2.9 to Flint 3.
- [gh-100](https://github.com/flintlib/python-flint/issues/100):
  Supports Python 3.12 by using setuptools instead of
  numpy.distutils.

New features:

- [gh-87](https://github.com/flintlib/python-flint/issues/87):
  Adds `fmpz_mod_poly` type for polynomials over `fmpz_mod`.
- [gh-85](https://github.com/flintlib/python-flint/issues/85):
  Adds discrete logarithms to `fmpz_mod`.
- [gh-83](https://github.com/flintlib/python-flint/issues/83):
  Introduces the `fmpz_mod` type for multi-precision integer mods.

Bug fixes:

- [gh-93](https://github.com/flintlib/python-flint/issues/93):
  Fixes a bug with `pow(int, int, fmpz)` which previously gave
  incorrect results.
- [gh-78](https://github.com/flintlib/python-flint/issues/78),
  [gh-79](https://github.com/flintlib/python-flint/issues/79):
  minor fixes for the `nmod` type.

0.4.4
-----

- [gh-75](https://github.com/flintlib/python-flint/issues/75),
  [gh-77](https://github.com/flintlib/python-flint/issues/77):
  finish bulk of the work in refactoring `python-flint` into
  submodules
- [gh-72](https://github.com/flintlib/python-flint/issues/72):
  The roots method of `arb_poly` is not supported. Use either the
  `complex_roots` method or `acb_roots(p).roots()` to get the old
  behaviour of returning the complex roots. The `roots` method on
  `fmpz_poly` and `fmpq_poly` now return integer and rational
  roots respectively. To access complex roots on these types, use
  the `complex_roots` method. For `acb_poly`, both `roots` and
  `complex_roots` behave the same
- [gh-71](https://github.com/flintlib/python-flint/issues/71):
  Include files in sdist and fix issue
  [gh-70](https://github.com/flintlib/python-flint/issues/70)
- [gh-67](https://github.com/flintlib/python-flint/issues/67):
  Continue refactoring job to introduce submodules into `python-flint`

0.4.3
-----

- [gh-63](https://github.com/flintlib/python-flint/issues/63):
  The `roots` method of `arb_poly`, and `nmod_poly` is no longer
  supported. Use `acb_roots(p).roots()` to get the old behaviour
  of returning the roots as `acb`. Note that the `roots` method of
  `fmpz_poly` and `fmpq_poly` currently returns the complex roots
  of the polynomial.
- [gh-61](https://github.com/flintlib/python-flint/issues/61):
  Start refactoring job to introduce submodules into `python-flint`

0.4.2
-----

- [gh-57](https://github.com/flintlib/python-flint/issues/57):
  Adds manylinux wheels

0.4.1
-----

- [gh-47](https://github.com/flintlib/python-flint/issues/47):
  Removes Linux wheels, updates instructions for building from
  source.

0.4.0
-----

- [gh-45](https://github.com/flintlib/python-flint/issues/45):
  Adds wheels for Windows, OSX and manylinux but the Linux wheels
  are broken.

License
------------

Python-FLINT is licensed MIT. FLINT and Arb are LGPL v2.1+.
