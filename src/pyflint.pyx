"""
Python wrapper for FLINT - Fast Library for Number Theory
http://www.flintlib.org/
"""

cimport flint
cimport libc.stdlib
cimport cython

cdef flint_rand_t global_random_state
flint_randinit(global_random_state)

cdef extern from "Python.h":
    ctypedef void PyObject
    ctypedef void PyTypeObject
    ctypedef long Py_ssize_t
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

def matrix_to_str(tab):
    if len(tab) == 0 or len(tab[0]) == 0:
        return "[]"
    tab = [map(str, row) for row in tab]
    widths = []
    for i in xrange(len(tab[0])):
        w = max([len(row[i]) for row in tab])
        widths.append(w)
    for i in xrange(len(tab)):
        tab[i] = [s.rjust(widths[j]) for j, s in enumerate(tab[i])]
        tab[i] = "[" + (", ".join(tab[i])) + "]"
    return "\n".join(tab)

cdef inline bint typecheck(object ob, object tp):
    return PyObject_TypeCheck(ob, <PyTypeObject*>tp)

DEF FMPZ_UNKNOWN = 0
DEF FMPZ_REF = 1
DEF FMPZ_TMP = 2


cdef long prec_to_dps(n):
    return max(1, int(round(int(n)/3.3219280948873626)-1))

cdef long dps_to_prec(n):
    return max(1, int(round((int(n)+1)*3.3219280948873626)))

cdef class Context:
    cpdef public bint pretty
    cpdef public long _prec
    cpdef public long _dps
    cpdef arf_rnd_t rnd

    def __init__(self):
        self.default()

    def default(self):
        self.pretty = False
        self.rnd = ARF_RND_DOWN
        self.prec = 53

    property prec:

        def __set__(self, prec):
            cdef long cprec = prec
            if cprec < 2:
                raise ValueError("prec must be >= 2")
            self._prec = cprec
            self._dps = prec_to_dps(cprec)

        def __get__(self):
            return self._prec

    property dps:

        def __set__(self, prec):
            self.prec = dps_to_prec(prec)

        def __get__(self):
            return self._dps


    def __repr__(self):
        return "pretty = %s\nprec = %s" % (self.pretty, self.prec)

cdef Context thectx = Context()

ctx = thectx

cdef inline long getprec(long prec=0):
    if prec > 1:
        return prec
    else:
        return thectx._prec

include "fmpz.pyx"
include "fmpz_poly.pyx"
include "fmpz_mat.pyx"

include "fmpq.pyx"
include "fmpq_poly.pyx"
include "fmpq_mat.pyx"

include "nmod.pyx"
include "nmod_poly.pyx"
include "nmod_mat.pyx"

include "arf.pyx"
include "arb.pyx"
include "arb_poly.pyx"
include "arb_mat.pyx"

include "acb.pyx"
include "acb_poly.pyx"
include "acb_mat.pyx"

include "functions.pyx"

