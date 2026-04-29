from flint.flint_base.flint_base cimport flint_poly

from flint.flintlib.functions.acb_poly cimport acb_poly_t

cdef acb_poly_set_list(acb_poly_t poly, list val, long prec)
cdef class acb_poly(flint_poly):
    cdef acb_poly_t val
    cpdef long length(self)
    cpdef long degree(self)
