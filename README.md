Python-FLINT
============

Python extension module wrapping FLINT (Fast Library for Number Theory)
and Arb (arbitrary-precision ball arithmetic). Features:

* Integers, rationals, integers mod n
* Real and complex numbers with rigorous error tracking
* Polynomials and matrices over all the above types
* Lots of special functions

Author: Fredrik Johansson <fredrik.johansson@gmail.com>

Repository: https://github.com/fredrik-johansson/python-flint/

Installation
------------

First install both FLINT and Arb (currently, the git versions are required).
See:

* https://github.com/fredrik-johansson/flint2/
* https://github.com/fredrik-johansson/arb/

Install build dependencies:

    sudo apt-get install cython python-dev

Build and install python-flint from source:

    python setup.py build_ext
    sudo python setup.py install

Run the test suite:

    python test/test.py

Import using:

    >>> from flint import *

### Additional paths

The FLINT and Arb header files and library files (libflint.so and libarb.so)
must be available at compile time. If they are in a nonstandard location
(for example, if they have been built but not installed),
use a command such as the following to build:

    python ./setup.py build_ext \
        --include-dirs=/home/fredrik/src/flint2:/home/fredrik/src/arb \
        --library-dirs=/home/fredrik/src/flint2:/home/fredrik/src/arb

Likewise, before starting the Python interpreter, tell the linker
where to find the library files using something like:

    export LD_LIBRARY_PATH=/home/fredrik/src/flint2:/home/fredrik/src/arb:$LD_LIBRARY_PATH

You may also have to install the CPimport file:

    sudo cp /home/fredrik/src/flint2/qadic/CPimport.txt /usr/local/share/flint/CPimport.txt

Examples
-------------------------------------

Number-theoretic functions:

    >>> fmpz(1000).number_of_partitions()
    fmpz(24061467864032622473692149727991)
    >>> fmpq.bernoulli_ui(64)
    fmpq(-106783830147866529886385444979142647942017,510)

Polynomial arithmetic:

    >>> a = fmpz_poly([1,2,3]); b = fmpz_poly([2,3,4]); a.gcd(a * b)
    fmpz_poly([1, 2, 3])
    >>> a = fmpz_poly(range(10001)); b = fmpz_poly(range(10000)); a.gcd(a * b).degree()
    10000

Matrix arithmetic:

    >>> (fmpz_mat(2,2,[1,1,1,0]) ** 10)
    fmpz_mat(2, 2, [89, 55, 55, 34])

Numerical evaluation:

    >>> showgood(lambda: (arb.pi() * arb(163).sqrt()).exp() - 640320**3 - 744, dps=25)
    -7.499274028018143111206461e-13
    >>> showgood(lambda: (arb.pi() * 10**100 + arb(1)/1000).sin(), dps=25)
    0.0009999998333333416666664683


Documentation
-------------------------------------

Complete API documentation is available.

To do
-------------------------------------

* Power series
* Finite fields
* p-adic numbers
* Rational functions
* Vector or array types (maybe)
* Many convenience methods
* Lots of FLINT and Arb functions/methods that haven't been wrapped yet
* Write generic implementations of functions missing for specific FLINT types
* Proper handling of special values in various places (throwing Python exceptions instead of aborting, etc.)
* Various automatic conversions
* Conversions to and from external types (numpy, sage, sympy, mpmath, gmpy)
* Support for Python-level multithreading (there is some global state)
* Support for subclassing (maybe, not a priority)
* Improved printing and string input/output
* IPython hooks (TeX pretty-printing etc.)
* Python 3.x support
* Windows support

License
------------

Python-FLINT is licensed MIT. FLINT and Arb are GPL v2+.

