from flint.flint_base.flint_base cimport flint_scalar
from flint.utils.conversion cimport chars_from_str

from flint._flint cimport fmpz_t, slong, pylong_as_slong
from flint._flint cimport PyObject, fmpz_set_str, fmpz_set_si
from flint._flint cimport PyInt_Check, PyInt_AS_LONG, PyLong_Check

from cpython.version cimport PY_MAJOR_VERSION

cdef int fmpz_set_any_ref(fmpz_t x, obj)

cdef inline int fmpz_set_pylong(fmpz_t x, obj):
    cdef int overflow
    cdef slong longval
    longval = pylong_as_slong(<PyObject*>obj, &overflow)
    if overflow:
        s = "%x" % obj
        fmpz_set_str(x, chars_from_str(s), 16)
    else:
        fmpz_set_si(x, longval)

cdef inline int fmpz_set_python(fmpz_t x, obj):
    if PY_MAJOR_VERSION < 3 and PyInt_Check(<PyObject*>obj):
        fmpz_set_si(x, PyInt_AS_LONG(<PyObject*>obj))
        return 1
    if PyLong_Check(<PyObject*>obj):
        fmpz_set_pylong(x, obj)
        return 1
    return 0
cdef any_as_fmpz(obj)

cdef class fmpz(flint_scalar):
    """
    The *fmpz* type represents an arbitrary-size integer.

        >>> fmpz(3) ** 25
        847288609443

    """ 

    cdef fmpz_t val
