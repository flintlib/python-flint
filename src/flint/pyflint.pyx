"""
Python wrapper for FLINT and Arb.
"""

cimport flint
cimport libc.stdlib
cimport cython

from flint._global_context cimport thectx

cdef flint_rand_t global_random_state
flint_randinit(global_random_state)

cdef extern from "Python.h":
    int PyObject_TypeCheck(object, PyTypeObject*)
    int PyInt_Check(PyObject *o)
    PyObject* PyInt_FromLong(long ival)
    int PyLong_Check(PyObject *o)
    long PyInt_AS_LONG(PyObject *io)
    double PyFloat_AS_DOUBLE(PyObject *io)
    Py_ssize_t PyList_GET_SIZE(PyObject *list)
    long PyLong_AsLongAndOverflow(PyObject *pylong, int *overflow)
    int PyComplex_Check(PyObject *o)
    double PyComplex_RealAsDouble(PyObject *op)
    double PyComplex_ImagAsDouble(PyObject *op)

cdef inline bint typecheck(object ob, object tp):
    return PyObject_TypeCheck(ob, <PyTypeObject*>tp)

DEF FMPZ_UNKNOWN = 0
DEF FMPZ_REF = 1
DEF FMPZ_TMP = 2

ctx = thectx

# TODO:
# This should be factored out into flint_base
# but we cannot do this until we can import 
# acb_poly to allow the roots() method to remain

from flint.flint_base.flint_base cimport flint_elem

cdef class flint_poly(flint_elem):
    """
    Base class for polynomials.
    """

    def __iter__(self):
        cdef long i, n
        n = self.length()
        for i in range(n):
            yield self[i]

    def coeffs(self):
        return list(self)

    def str(self, bint ascending=False):
        """
        Convert to a human-readable string (generic implementation for
        all polynomial types).

        If *ascending* is *True*, the monomials are output from low degree to
        high, otherwise from high to low.
        """
        coeffs = [str(c) for c in self]
        if not coeffs:
            return "0"
        s = []
        coeffs = enumerate(coeffs)
        if not ascending:
            coeffs = reversed(list(coeffs))
        for i, c in coeffs:
            if c == "0":
                continue
            else:
                if c.startswith("-") or (" " in c):
                    c = "(" + c + ")"
                if i == 0:
                    s.append("%s" % c)
                elif i == 1:
                    if c == "1":
                        s.append("x")
                    else:
                        s.append("%s*x" % c)
                else:
                    if c == "1":
                        s.append("x^%s" % i)
                    else:
                        s.append("%s*x^%s" % (c, i))
        return " + ".join(s)

    def roots(self, **kwargs):
        """
        Isolates the complex roots of *self*. See :meth:`.acb_poly.roots`
        for details.
        """
        return acb_poly(self).roots(**kwargs)

include "fmpz.pyx"
include "fmpz_poly.pyx"
include "fmpz_mpoly.pyx"
include "fmpz_mat.pyx"
include "fmpz_series.pyx"

include "fmpq.pyx"
include "fmpq_poly.pyx"
include "fmpq_mat.pyx"
include "fmpq_series.pyx"

include "nmod.pyx"
include "nmod_poly.pyx"
include "nmod_mat.pyx"
include "nmod_series.pyx"

include "arf.pyx"
include "arb.pyx"
include "arb_poly.pyx"
include "arb_mat.pyx"
include "arb_series.pyx"

include "acb.pyx"
include "acb_poly.pyx"
include "acb_mat.pyx"
include "acb_series.pyx"

include "functions.pyx"

include "dirichlet.pyx"

