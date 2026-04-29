from flint.flint_base.flint_base cimport flint_mpoly
from flint.flint_base.flint_base cimport flint_mpoly_context

from flint.flintlib.functions.fmpq_mpoly cimport (
    fmpq_mpoly_ctx_t,
    fmpq_mpoly_t,
    fmpq_mpoly_struct,
    fmpq_mpoly_init,
)
from flint.flintlib.types.flint cimport slong

from flint.types.fmpz_mpoly cimport fmpz_mpoly_ctx

cdef inline init_fmpq_mpoly(fmpq_mpoly var, fmpq_mpoly_ctx ctx):
    var.ctx = ctx
    fmpq_mpoly_init(var.val, ctx.val)
    var._init = True

cdef inline create_fmpq_mpoly(fmpq_mpoly_ctx ctx):
    cdef fmpq_mpoly var
    var = fmpq_mpoly.__new__(fmpq_mpoly)
    var.ctx = ctx
    fmpq_mpoly_init(var.val, ctx.val)
    var._init = True
    return var

cdef class fmpq_mpoly_ctx(flint_mpoly_context):
    cdef fmpq_mpoly_ctx_t val

cdef class fmpq_mpoly(flint_mpoly):
    cdef fmpq_mpoly_t val
    cdef fmpq_mpoly_ctx ctx
    cdef bint _init

cdef class fmpq_mpoly_vec:
    cdef fmpq_mpoly_struct *val
    cdef slong length
    cdef fmpq_mpoly_ctx ctx
    cdef fmpq_mpoly_struct **double_indirect
