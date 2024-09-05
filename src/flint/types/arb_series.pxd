from flint.flint_base.flint_base cimport flint_series

from flint.flintlib.functions.arb_poly cimport arb_poly_t

cdef class arb_series(flint_series):
    cdef arb_poly_t val
    cdef long prec
    cpdef long length(self)
    cpdef valuation(self)
