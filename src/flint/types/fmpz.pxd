from cpython.long cimport PyLong_Check
from flint.flint_base.flint_base cimport flint_scalar
from flint.utils.conversion cimport chars_from_str
from flint.flintlib.types.flint cimport slong, ulong, pylong_as_slong
from flint.flintlib.types.flint cimport PyObject
from flint.flintlib.functions.fmpz cimport fmpz_t, fmpz_set_si, fmpz_set_signed_ui_array
import sys

cdef int fmpz_set_any_ref(fmpz_t x, obj)
cdef fmpz_get_intlong(fmpz_t x)

cdef int is_big_endian = int(sys.byteorder == "big")

cdef inline ulong ulong_from_little_endian(unsigned char *ptr):
    # Read a ulong from little-endian bytes
    cdef ulong w = 0
    for i in range(sizeof(ulong) // 8):
        w = (w << 8) | ptr[i]
    return w

cdef inline int fmpz_set_pylong(fmpz_t x, obj):
    cdef int overflow
    cdef slong longval
    cdef slong size
    cdef bytes b
    cdef ulong w
    cdef ulong *words
    cdef int i

    longval = pylong_as_slong(<PyObject*>obj, &overflow)
    if overflow:
        # make sure the sign bit fits
        # we need 8 * sizeof(ulong) * size > obj.bit_length()
        size = obj.bit_length() // (8 * sizeof(ulong)) + 1
        b = obj.to_bytes(sizeof(ulong) * size, "little", signed=True)
        # b is a local Python object, we access the internal pointer
        words = <ulong*>(<char *>b)
        if is_big_endian:
            for i in range(size):
                words[i] = ulong_from_little_endian(<unsigned char *>(words + i))
        fmpz_set_signed_ui_array(x, words, size)
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
