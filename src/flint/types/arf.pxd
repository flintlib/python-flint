from flint.flintlib.functions.arf cimport arf_t

cdef class arf:
    cdef arf_t val
    cpdef bint is_finite(self)
    cpdef bint is_pos_inf(self)
    cpdef bint is_neg_inf(self)
    cpdef bint is_nan(self)
    cpdef bint is_zero(self)
