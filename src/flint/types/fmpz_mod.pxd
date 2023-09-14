from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.fmpz_mod cimport fmpz_mod_ctx_t


cdef class fmpz_mod_ctx:
    cdef fmpz_mod_ctx_t val

cdef class fmpz_mod(flint_scalar):
    cdef fmpz_mod_ctx ctx
    cdef fmpz_t val

