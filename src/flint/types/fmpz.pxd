from cpython.long cimport PyLong_Check
from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.types.flint cimport PyObject, PyLongLayout, PyLong_GetNativeLayout, PyLongExport, PyLong_Export, PyLong_FreeExport
from flint.flintlib.functions.fmpz cimport fmpz_t, fmpz_set_str, fmpz_set_si, _fmpz_promote_val, _fmpz_demote_val, mpz_import, fmpz_neg, fmpz_ui_pow_ui, fmpz_init, fmpz_sub, fmpz_clear

cdef extern from "<limits.h>":
    const long LONG_MAX
    const long LONG_MIN

cdef int fmpz_set_any_ref(fmpz_t x, obj)
cdef fmpz_get_intlong(fmpz_t x)

cdef inline int fmpz_set_pylong(fmpz_t x, obj):
    cdef const PyLongLayout * layout
    cdef PyLongExport long_export
    cdef fmpz_t tmp

    layout = PyLong_GetNativeLayout()

    PyLong_Export(obj, &long_export)
    if long_export.digits:
        z = _fmpz_promote_val(x)
        mpz_import(z, long_export.ndigits, layout[0].digits_order,
                   layout[0].digit_size, layout[0].digit_endianness,
                   layout[0].digit_size*8 - layout[0].bits_per_digit,
                   long_export.digits)
        _fmpz_demote_val(x)
        if long_export.negative:
            fmpz_neg(x, x)
        PyLong_FreeExport(&long_export)
    else:
        if LONG_MIN <= long_export.value <= LONG_MAX:
            fmpz_set_si(x, long_export.value)
        else:
            z = _fmpz_promote_val(x)
            mpz_import(z, 1, -1, 8, 0, 0, &long_export.value);
            _fmpz_demote_val(x)
            if long_export.value < 0:
                fmpz_init(tmp)
                fmpz_ui_pow_ui(tmp, 2, 64)
                fmpz_sub(x, x, tmp)
                fmpz_clear(tmp)

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
