from flint.flint_base.flint_base cimport flint_mpoly, flint_mod_mpoly_context

from flint.flintlib.functions.fmpz_mod_mpoly cimport (
    fmpz_mod_mpoly_ctx_t,
    fmpz_mod_mpoly_t,
    fmpz_mod_mpoly_init,
    fmpz_mod_mpoly_struct
)
from flint.flintlib.types.flint cimport slong

cdef inline init_fmpz_mod_mpoly(fmpz_mod_mpoly var, fmpz_mod_mpoly_ctx ctx):
    var.ctx = ctx
    fmpz_mod_mpoly_init(var.val, ctx.val)
    var._init = True

cdef inline fmpz_mod_mpoly create_fmpz_mod_mpoly(fmpz_mod_mpoly_ctx ctx):
    cdef fmpz_mod_mpoly var
    var = fmpz_mod_mpoly.__new__(fmpz_mod_mpoly)
    var.ctx = ctx
    fmpz_mod_mpoly_init(var.val, ctx.val)
    var._init = True
    return var

cdef class fmpz_mod_mpoly_ctx(flint_mod_mpoly_context):
    cdef fmpz_mod_mpoly_ctx_t val

cdef class fmpz_mod_mpoly(flint_mpoly):
    cdef fmpz_mod_mpoly_t val
    cdef fmpz_mod_mpoly_ctx ctx
    cdef bint _init

cdef class fmpz_mod_mpoly_vec:
    cdef fmpz_mod_mpoly_struct *val
    cdef slong length
    cdef fmpz_mod_mpoly_ctx ctx
    cdef fmpz_mod_mpoly_struct **double_indirect
