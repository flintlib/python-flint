"""
Python wrapper for FLINT and Arb.
"""

cimport flint
cimport libc.stdlib
cimport cython

from flint.flint_base.flint_context cimport thectx

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

#DEF FMPZ_UNKNOWN = 0
#DEF FMPZ_REF = 1
#DEF FMPZ_TMP = 2

ctx = thectx

# include "fmpz.pyx"
# include "fmpz_poly.pyx"
include "fmpz_mpoly.pyx"
# include "fmpz_mat.pyx"
# include "fmpz_series.pyx"

# include "fmpq.pyx"
# include "fmpq_poly.pyx"
# include "fmpq_mat.pyx"
# include "fmpq_series.pyx"

# include "nmod.pyx"
# include "nmod_poly.pyx"
# include "nmod_mat.pyx"
# include "nmod_series.pyx"

# include "arf.pyx"
# include "arb.pyx"
# include "arb_poly.pyx"
# include "arb_mat.pyx"
# include "arb_series.pyx"

# include "acb.pyx"
include "acb_poly.pyx"
include "acb_mat.pyx"
include "acb_series.pyx"

include "functions.pyx"

include "dirichlet.pyx"

