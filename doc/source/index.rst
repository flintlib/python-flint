.. python-flint documentation master file, created by
   sphinx-quickstart on Wed Jul 18 11:49:46 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Python-FLINT
==========================================================

Python extension module wrapping FLINT_ (Fast Library for Number Theory)
and Arb_ (arbitrary-precision ball arithmetic). Features:

.. _FLINT: http://flintlib.org/
.. _Arb: http://arblib.org/


* Integers, rationals, integers mod n
* Real and complex numbers with rigorous error tracking
* Polynomials and matrices over all the above types
* Lots of mathematical functions

Author: `Fredrik Johansson <http://fredrikj.net/>`_ <fredrik.johansson@gmail.com>

Repository: https://github.com/fredrik-johansson/python-flint/

Introduction
------------

.. toctree::
   :maxdepth: 2

   setup.rst
   general.rst


Reference
---------

Scalar types
............

.. toctree::
   :maxdepth: 1

   fmpz.rst
   fmpq.rst
   nmod.rst
   arb.rst
   acb.rst
   dirichlet.rst

Matrix types
............

.. toctree::
   :maxdepth: 1

   fmpz_mat.rst
   fmpq_mat.rst
   nmod_mat.rst
   arb_mat.rst
   acb_mat.rst

Polynomial types
................

.. toctree::
   :maxdepth: 1

   fmpz_poly.rst
   fmpq_poly.rst
   nmod_poly.rst
   arb_poly.rst
   acb_poly.rst

Power series types
..................

.. toctree::
   :maxdepth: 1

   fmpz_series.rst
   fmpq_series.rst
   nmod_series.rst
   arb_series.rst
   acb_series.rst

