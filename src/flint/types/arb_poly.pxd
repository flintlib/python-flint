from flint.flint_base.flint_base cimport flint_poly
from flint.flintlib.functions.arb_poly cimport arb_poly_t

cdef arb_poly_set_list(arb_poly_t poly, list val, long prec)

cdef class arb_poly(flint_poly):
    cdef arb_poly_t val
    cpdef long length(self)
    cpdef long degree(self)
