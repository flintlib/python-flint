from flint.flint_base.flint_base cimport flint_mpoly

from flint.flintlib.fmpz_mpoly cimport fmpz_mpoly_ctx_t
from flint.flintlib.fmpz_mpoly cimport fmpz_mpoly_t
from flint.flintlib.flint cimport slong

cdef class fmpz_mpoly_ctx:
    cdef fmpz_mpoly_ctx_t val
    cpdef slong nvars(self)

    cpdef ordering(self)

cdef class fmpz_mpoly(flint_mpoly):
    cdef fmpz_mpoly_t val
    cdef fmpz_mpoly_ctx ctx
    cdef bint _init
