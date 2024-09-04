from flint.flint_base.flint_base cimport flint_series

from flint.flintlib.functions.fmpz_poly cimport fmpz_poly_t

cdef class fmpz_series(flint_series):
    cdef fmpz_poly_t val
    cdef long prec
    cpdef long length(self)
    cpdef valuation(self)
