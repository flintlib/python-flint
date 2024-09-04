from flint.flint_base.flint_base cimport flint_poly
from flint.flintlib.functions.fmpq_poly cimport fmpq_poly_t

cdef fmpq_poly_set_list(fmpq_poly_t poly, list val)
cdef any_as_fmpq_poly(obj)
cdef class fmpq_poly(flint_poly):
    cdef fmpq_poly_t val
    cpdef long length(self)
    cpdef long degree(self)
