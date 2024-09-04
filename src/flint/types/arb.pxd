from flint.flint_base.flint_base cimport flint_scalar

from flint.flintlib.functions.arb cimport arb_t

cdef any_as_arb_or_notimplemented(x)
cdef int arb_set_python(arb_t x, obj, bint allow_conversion) except -1
cdef any_as_arb(x)
cdef arb_set_mpmath_mpf(arb_t x, obj)

cdef class arb(flint_scalar):
    cdef arb_t val

    cpdef bint is_zero(self)
    cpdef bint is_finite(self)
    cpdef bint is_nan(self)
    cpdef bint is_exact(self)
    cpdef bint is_integer(self)
