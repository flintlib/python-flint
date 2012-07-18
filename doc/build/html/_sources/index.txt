.. python-flint documentation master file, created by
   sphinx-quickstart on Wed Jul 18 11:49:46 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Python bindings for FLINT
==========================================================

Introduction
::::::::::::

The python-flint Python extension module provides
bindings to FLINT (http://flintlib.org/).
FLINT (Fast Library for Number Theory) is a C library that implements
highly optimised exact scalar, polynomial and matrix arithmetic as
well as many number-theoretic special functions.

As of the current version, only a subset of the FLINT types
are wrapped, and only a very limited selection of methods are provided.
There are known and unknown bugs, and the interface is subject to
change drastically in the future. Any feedback is welcome
and can be sent to fredrik.johansson@gmail.com or flint-devel@googlegroups.com.
Of particular note:

* Most functionality in FLINT is not yet wrapped!
* Coercions are not supported in all directions.
* Conversion between Python longs and fmpzs is not done efficiently.
* Types don't yet provide hash methods.
* Interaction with gmpy (http://gmpy.org/) types (and types from many other libraries) is not supported.

Setup
::::::::::::

To build this module, you first need to install the latest
version of FLINT (follow the instructions in the FLINT documentation). You
also need to install Cython (http://cython.org/).

Next, extract the python-flint archive and run::

    python ./setup.py build_ext
    sudo python ./setup.py install
    python ./test.py

Note: if FLINT is not installed globally, you can do something like the following (the second
line is required for loading FLINT from Python)::

    python ./setup.py build_ext --include-dirs=/home/fredrik/src/flint2 --library-dirs=/home/fredrik/src/flint2
    export LD_LIBRARY_PATH=/home/fredrik/src/flint2/:$LD_LIBRARY_PATH

If everything worked, you will be able to import the ``flint`` module
in Python::

    Python 2.7.3 (default, Apr 20 2012, 22:39:59) 
    [GCC 4.6.3] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import flint
    >>> flint.fmpz(33) ** 33
    fmpz(129110040087761027839616029934664535539337183380513)

Reference
::::::::::::

Scalar types
............

fmpz
------------------------------------------------------------

.. autoclass :: flint.fmpz
  :members:
  :undoc-members:

fmpq
------------------------------------------------------------

.. autoclass :: flint.fmpq
  :members:
  :undoc-members:

nmod
------------------------------------------------------------

.. autoclass :: flint.nmod
  :members:
  :undoc-members:

Polynomial types
................

fmpz_poly
------------------------------------------------------------

.. autoclass :: flint.fmpz_poly
  :members:
  :undoc-members:

fmpq_poly
------------------------------------------------------------

.. autoclass :: flint.fmpq_poly
  :members:
  :undoc-members:

nmod_poly
------------------------------------------------------------

.. autoclass :: flint.nmod_poly
  :members:
  :undoc-members:

Matrix types
............

fmpz_mat
------------------------------------------------------------

.. autoclass :: flint.fmpz_mat
  :members:
  :undoc-members:

fmpq_mat
------------------------------------------------------------

.. autoclass :: flint.fmpq_mat
  :members:
  :undoc-members:

nmod_mat
------------------------------------------------------------

.. autoclass :: flint.nmod_mat
  :members:
  :undoc-members:

Special functions
.................

.. autofunction :: flint.number_of_partitions
