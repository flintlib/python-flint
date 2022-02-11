"""
Python wrapper for FLINT and Arb.
"""

cimport flint
cimport libc.stdlib
cimport cython

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

from cpython.version cimport PY_MAJOR_VERSION

cdef chars_from_str(s):
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.encode('ascii')

cdef str_from_chars(s):
    if PY_MAJOR_VERSION < 3:
        return str(s)
    else:
        return bytes(s).decode('ascii')

cdef matrix_to_str(tab):
    if len(tab) == 0 or len(tab[0]) == 0:
        return "[]"
    tab = [[str(c) for c in row] for row in tab]
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

cdef class FlintContext:
    cdef public bint pretty
    cdef public long _prec
    cdef public long _dps
    cdef arf_rnd_t rnd
    cdef public bint unicode
    cdef public long _cap

    def __init__(self):
        self.default()

    def default(self):
        self.pretty = True
        self.rnd = ARF_RND_DOWN
        self.prec = 53
        self.unicode = False
        self.threads = 1
        self.cap = 10

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

    property cap:

        def __set__(self, long cap):
            if cap < 0:
                raise ValueError("cap must be >= 0")
            self._cap = cap

        def __get__(self):
            return self._cap

    property threads:

        def __set__(self, long num):
            assert num >= 1 and num <= 64
            flint_set_num_threads(num)

        def __get__(self):
            return flint_get_num_threads()

    def __repr__(self):
        return "pretty = %-8s  # pretty-print repr() output\n" \
               "unicode = %-8s # use unicode characters in output\n" \
               "prec = %-8s    # real/complex precision (in bits)\n"   \
               "dps = %-8s     # real/complex precision (in digits)\n"    \
               "cap = %-8s     # power series precision\n"    \
               "threads = %-8s # max number of threads used internally\n" % \
            (self.pretty, self.unicode, self.prec, self.dps, self.cap, self.threads)

    def cleanup(self):
        flint_cleanup()

cdef FlintContext thectx = FlintContext()

ctx = thectx

cdef inline long getprec(long prec=0):
    if prec > 1:
        return prec
    else:
        return thectx._prec

cdef inline long getcap(long cap=-1):
    if cap >= 0:
        return cap
    else:
        return thectx._cap

cdef class flint_elem:

    def __repr__(self):
        if ctx.pretty:
            return self.str()
        else:
            return self.repr()

    def __str__(self):
        return self.str()

cdef class flint_scalar(flint_elem):
    pass

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

cdef class flint_mpoly(flint_elem):
    """
    Base class for multivariate polynomials.
    """


cdef class flint_mat(flint_elem):
    """
    Base class for matrices.
    """

    def repr(self):
        if ctx.pretty:
            return str(self)
        # XXX
        return "%s(%i, %i, [%s])" % (type(self).__name__,
            self.nrows(), self.ncols(), (", ".join(map(str, self.entries()))))

    def str(self, *args, **kwargs):
        tab = self.table()
        if len(tab) == 0 or len(tab[0]) == 0:
            return "[]"
        tab = [[r.str(*args, **kwargs) for r in row] for row in tab]
        widths = []
        for i in xrange(len(tab[0])):
            w = max([len(row[i]) for row in tab])
            widths.append(w)
        for i in xrange(len(tab)):
            tab[i] = [s.rjust(widths[j]) for j, s in enumerate(tab[i])]
            tab[i] = "[" + (", ".join(tab[i])) + "]"
        return "\n".join(tab)

    def entries(self):
        cdef long i, j, m, n
        m = self.nrows()
        n = self.ncols()
        L = [None] * (m * n)
        for i from 0 <= i < m:
            for j from 0 <= j < n:
                L[i*n + j] = self[i, j]
        return L

    def __iter__(self):
        cdef long i, j, m, n
        m = self.nrows()
        n = self.ncols()
        for i from 0 <= i < m:
            for j from 0 <= j < n:
                yield self[i, j]

    def table(self):
        cdef long i, m, n
        m = self.nrows()
        n = self.ncols()
        L = self.entries()
        return [L[i*n : (i+1)*n] for i in range(m)]

    # supports mpmath conversions
    tolist = table

cdef class flint_series(flint_elem):
    """
    Base class for power series.
    """
    def __iter__(self):
        cdef long i, n
        n = self.length()
        for i in range(n):
            yield self[i]

    def coeffs(self):
        return list(self)


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

