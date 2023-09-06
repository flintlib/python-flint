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

For Windows (x86-64) or OSX (x86-64 or arm64) or Linux (x86-64
`manylinux_2_17`) there are binary wheels for python-flint on PyPI. For these
platforms python-flint can be installed simply with `pip`

    pip install python-flint

Alternatively python-flint can be installed using `conda`

    conda install -c conda-forge python-flint

Build from source
-----------------

For other platforms or architectures installation needs to build from source.
First install both FLINT 2.9.0 and Arb 2.23. Note that for python-flint 0.4
only these *exact* versions of FLINT and Arb will work. While some Linux
distributions may provide FLINT and Arb it is unlikely that they will provide
the exact versions required (e.g. for Ubuntu only Ubuntu 23.04 provides these
versions at the time of writing).

See here for instructions on building FLINT and Arb:

* http://flintlib.org/
* http://arblib.org/

The latest release of Python-FLINT can then be built and installed using:

    pip install 'cython>=3' numpy wheel
    pip install --no-build-isolation python-flint

Python-FLINT can also be installed from a git checkout or a source archive
as follows:

    pip install 'cython>=3' numpy wheel
    pip install --no-build-isolation .

A script that builds and installs FLINT, Arb and Python-FLINT that is tested on
Ubuntu can be found in the git repo here:

* https://github.com/flintlib/python-flint/blob/master/bin/pip_install_ubuntu.sh

See the documentation for further notes on building and installing
Python-FLINT:

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
* Wrap missing flint types: finite fields, p-adic numbers, multiprecision integer mods, rational functions
* Vector or array types (maybe)
* Many convenience methods
* Write generic implementations of functions missing for specific FLINT types
* Proper handling of special values in various places (throwing Python exceptions instead of aborting, etc.)
* Various automatic conversions
* Conversions to and from external types (numpy, sage, sympy, mpmath, gmpy)
* Improved printing and string input/output
* IPython hooks (TeX pretty-printing etc.)

License
------------

Python-FLINT is licensed MIT. FLINT and Arb are LGPL v2.1+.
