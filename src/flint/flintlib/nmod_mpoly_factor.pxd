from flint.flintlib.nmod_mpoly cimport nmod_mpoly_ctx_t, nmod_mpoly_struct, nmod_mpoly_t
from flint.flintlib.types.flint cimport slong, ulong, mp_limb_t
from flint.flintlib.fmpz cimport fmpz_struct


cdef extern from "flint/nmod_mpoly_factor.h":
    ctypedef struct nmod_mpoly_factor_struct:
        mp_limb_t constant
        nmod_mpoly_struct * poly
        fmpz_struct * exp
        slong num
        slong alloc

    ctypedef nmod_mpoly_factor_struct nmod_mpoly_factor_t[1]

    void nmod_mpoly_factor_init(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_factor_clear(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_factor_swap(nmod_mpoly_factor_t f, nmod_mpoly_factor_t g, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_factor_length(const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    ulong nmod_mpoly_factor_get_constant_ui(const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_factor_get_base(nmod_mpoly_t p, const nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_factor_swap_base(nmod_mpoly_t p, nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_factor_get_exp_si(nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_factor_sort(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_factor_squarefree(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_factor(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
