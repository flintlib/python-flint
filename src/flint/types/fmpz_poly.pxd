from flint.flint_base.flint_base cimport flint_poly

from flint.flintlib.functions.fmpz_poly cimport fmpz_poly_t

cdef fmpz_poly_set_list(fmpz_poly_t poly, list val)

cdef any_as_fmpz_poly(x)

cdef class fmpz_poly(flint_poly):
    cdef fmpz_poly_t val
    cpdef long length(self)
    cpdef long degree(self)
