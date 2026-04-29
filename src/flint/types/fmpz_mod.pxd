from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.functions.fmpz cimport fmpz_struct, fmpz_t
from flint.flintlib.functions.fmpz_mod cimport (
    fmpz_mod_ctx_t,
    fmpz_mod_discrete_log_pohlig_hellman_t
)

cdef class fmpz_mod_ctx:
    cdef fmpz_mod_ctx_t val
    cdef bint _is_prime
    cdef bint _init_L
    cdef fmpz_mod_discrete_log_pohlig_hellman_t L

    cdef set_any_as_fmpz_mod(self, fmpz_t val, obj)
    cdef any_as_fmpz_mod(self, obj)
    cdef discrete_log_pohlig_hellman_run(self, fmpz_t x, fmpz_t y)

cdef class fmpz_mod(flint_scalar):
    cdef fmpz_mod_ctx ctx
    cdef fmpz_t val
    cdef fmpz_t x_g
