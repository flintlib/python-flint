from flint.flint_base.flint_base cimport flint_series
from flint.flintlib.functions.fmpq_poly cimport fmpq_poly_t

cdef class fmpq_series(flint_series):
    cdef fmpq_poly_t val
    cdef long prec
    cpdef long length(self)
    cpdef valuation(self)
    cdef bint zero_constant_term(s)
    cdef bint one_constant_term(s)
