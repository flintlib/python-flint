from flint.flint_base.flint_base cimport flint_mpoly
from flint.flint_base.flint_base cimport flint_mpoly_context

from flint.flintlib.fmpz_mpoly cimport fmpz_mpoly_ctx_t
from flint.flintlib.fmpz_mpoly cimport fmpz_mpoly_t, fmpz_mpoly_init
from flint.flintlib.flint cimport slong

cdef inline init_fmpz_mpoly(fmpz_mpoly var, fmpz_mpoly_ctx ctx):
    var.ctx = ctx
    fmpz_mpoly_init(var.val, ctx.val)
    var._init = True

cdef inline create_fmpz_mpoly(fmpz_mpoly_ctx ctx):
    cdef fmpz_mpoly var
    var = fmpz_mpoly.__new__(fmpz_mpoly)
    var.ctx = ctx
    fmpz_mpoly_init(var.val, ctx.val)
    var._init = True
    return var

cdef class fmpz_mpoly_ctx(flint_mpoly_context):
    cdef fmpz_mpoly_ctx_t val
    cpdef slong nvars(self)
    cpdef ordering(self)

cdef class fmpz_mpoly(flint_mpoly):
    cdef fmpz_mpoly_t val
    cdef fmpz_mpoly_ctx ctx
    cdef bint _init
