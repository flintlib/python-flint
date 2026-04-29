from flint.flintlib.functions.fmpz_mod_poly cimport *

from flint.flint_base.flint_base cimport flint_poly
from flint.types.fmpz_mod cimport fmpz_mod_ctx


cdef class fmpz_mod_poly_ctx:
    cdef fmpz_mod_ctx mod
    cdef any_as_fmpz_mod_poly(self, obj)
    cdef set_any_as_fmpz_mod_poly(self, fmpz_mod_poly_t poly, obj)
    cdef set_list_as_fmpz_mod_poly(self, fmpz_mod_poly_t poly, val)
    cdef new_ctype_poly(self)

cdef class fmpz_mod_poly(flint_poly):
    cdef fmpz_mod_poly_t val
    cdef fmpz_mod_poly_ctx ctx
    cpdef long length(self)
    cpdef long degree(self)
