from cpython.int cimport PyInt_Check, PyInt_AS_LONG
from cpython.long cimport PyLong_Check
from flint.flint_base.flint_base cimport flint_scalar
from flint.utils.conversion cimport chars_from_str

from flint.flintlib.flint cimport slong, pylong_as_slong
from flint.flintlib.flint cimport PyObject
from flint.flintlib.fmpz cimport fmpz_t, fmpz_set_str, fmpz_set_si

from cpython.version cimport PY_MAJOR_VERSION

cdef int fmpz_set_any_ref(fmpz_t x, obj)
cdef fmpz_get_intlong(fmpz_t x)

cdef int fmpz_set_pylong(fmpz_t x, obj)

cdef int fmpz_set_python(fmpz_t x, obj)

cdef any_as_fmpz(obj)

cdef class fmpz(flint_scalar):
    """
    The *fmpz* type represents an arbitrary-size integer.

        >>> fmpz(3) ** 25
        847288609443

    """ 

    cdef fmpz_t val
