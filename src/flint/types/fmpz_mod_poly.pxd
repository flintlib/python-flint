from flint.flintlib.fmpz_mod_poly cimport *

from flint.flint_base.flint_base cimport flint_poly
from flint.types.fmpz_mod cimport fmpz_mod_ctx

cdef class fmpz_mod_poly_ctx:
    cdef fmpz_mod_ctx mod

cdef class fmpz_mod_poly(flint_poly):
    cdef fmpz_mod_poly_t val
    cdef fmpz_mod_ctx ctx
    cpdef long length(self)
    cpdef long degree(self)
