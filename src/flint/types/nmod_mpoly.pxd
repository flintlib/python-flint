from flint.flint_base.flint_base cimport flint_mpoly, flint_mod_mpoly_context

from flint.flintlib.functions.nmod_mpoly cimport (
    nmod_mpoly_ctx_t,
    nmod_mpoly_t,
    nmod_mpoly_init,
    nmod_mpoly_struct
)
from flint.flintlib.types.flint cimport slong, ulong

cdef inline init_nmod_mpoly(nmod_mpoly var, nmod_mpoly_ctx ctx):
    var.ctx = ctx
    nmod_mpoly_init(var.val, ctx.val)
    var._init = True

cdef inline nmod_mpoly create_nmod_mpoly(nmod_mpoly_ctx ctx):
    cdef nmod_mpoly var
    var = nmod_mpoly.__new__(nmod_mpoly)
    var.ctx = ctx
    nmod_mpoly_init(var.val, ctx.val)
    var._init = True
    return var

cdef class nmod_mpoly_ctx(flint_mod_mpoly_context):
    cdef nmod_mpoly_ctx_t val

cdef class nmod_mpoly(flint_mpoly):
    cdef nmod_mpoly_t val
    cdef nmod_mpoly_ctx ctx
    cdef bint _init

cdef class nmod_mpoly_vec:
    cdef nmod_mpoly_struct *val
    cdef slong length
    cdef nmod_mpoly_ctx ctx
    cdef nmod_mpoly_struct **double_indirect
