from flint.flint_base.flint_base cimport flint_mpoly, flint_mpoly_context

from flint.flintlib.types.fmpz cimport (
    fmpz_mpoly_ctx_t,
    fmpz_mpoly_vec_t,
    fmpz_mpoly_t,
    fmpz_mpoly_struct,
)
from flint.flintlib.functions.fmpz_mpoly cimport fmpz_mpoly_init


cdef inline init_fmpz_mpoly(fmpz_mpoly var, fmpz_mpoly_ctx ctx):
    var.ctx = ctx
    fmpz_mpoly_init(var.val, ctx.val)
    var._init = True

cdef inline fmpz_mpoly create_fmpz_mpoly(fmpz_mpoly_ctx ctx):
    cdef fmpz_mpoly var
    var = fmpz_mpoly.__new__(fmpz_mpoly)
    var.ctx = ctx
    fmpz_mpoly_init(var.val, ctx.val)
    var._init = True
    return var

cdef class fmpz_mpoly_ctx(flint_mpoly_context):
    cdef fmpz_mpoly_ctx_t val

cdef class fmpz_mpoly(flint_mpoly):
    cdef fmpz_mpoly_t val
    cdef fmpz_mpoly_ctx ctx
    cdef bint _init

cdef class fmpz_mpoly_vec:
    cdef fmpz_mpoly_vec_t val
    cdef fmpz_mpoly_ctx ctx
    cdef fmpz_mpoly_struct **double_indirect
