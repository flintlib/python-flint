from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.functions.acb cimport acb_t

cdef any_as_acb(x)
cdef any_as_acb_or_notimplemented(x)
cdef int acb_set_python(acb_t x, obj, bint allow_conversion)
cdef class acb(flint_scalar):
    cdef acb_t val
    cpdef bint is_zero(self)
    cpdef bint is_finite(self)
    cpdef bint is_exact(self)
