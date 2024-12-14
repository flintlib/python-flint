from flint.flint_base.flint_base cimport flint_series

from flint.flintlib.functions.acb_poly cimport acb_poly_t

cdef class acb_series(flint_series):
    cdef acb_poly_t val
    cdef long prec
    cpdef long length(self)
    cpdef valuation(self)
