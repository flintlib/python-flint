from flint.flint_base.flint_base cimport flint_rational_function
from flint.flintlib.fmpz_mpoly_q cimport fmpz_mpoly_q_t

from flint.types.fmpz_mpoly cimport fmpz_mpoly, fmpz_mpoly_ctx

cdef class fmpz_mpoly_q(flint_rational_function):
    cdef fmpz_mpoly_q_t fraction
    cdef fmpz_mpoly_ctx ctx
    cdef bint _init
