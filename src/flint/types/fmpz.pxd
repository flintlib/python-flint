from cpython.long cimport PyLong_Check
from flint.flint_base.flint_base cimport flint_scalar
from flint.utils.conversion cimport chars_from_str
from flint.flintlib.types.flint cimport slong, pylong_as_slong
from flint.flintlib.types.flint cimport PyObject
from flint.flintlib.functions.fmpz cimport fmpz_t, fmpz_set_str, fmpz_set_si

cdef int fmpz_set_any_ref(fmpz_t x, obj)
cdef fmpz_get_intlong(fmpz_t x)

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
    if PyLong_Check(obj):
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
