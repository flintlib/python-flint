from flint.flint_base.flint_base cimport flint_mpoly
from flint.flint_base.flint_base cimport flint_mpoly_context

from flint.flintlib.fmpq_mpoly cimport fmpq_mpoly_ctx_t
from flint.flintlib.fmpq_mpoly cimport fmpq_mpoly_t
from flint.flintlib.flint cimport slong

cdef class fmpq_mpoly_ctx(flint_mpoly_context):
    cdef fmpq_mpoly_ctx_t val
    cpdef slong nvars(self)
    cpdef ordering(self)

cdef class fmpq_mpoly(flint_mpoly):
    cdef fmpq_mpoly_t val
    cdef fmpq_mpoly_ctx ctx
    cdef bint _init
