from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.fmpq cimport fmpq_mpoly_ctx_t, fmpq_mpoly_factor_t, fmpq_mpoly_t, fmpq_t



cdef extern from "flint/fmpq_mpoly_factor.h":
    void fmpq_mpoly_factor_init(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_clear(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    slong fmpq_mpoly_factor_length(const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_get_base(fmpq_mpoly_t B, const fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_swap_base(fmpq_mpoly_t B, fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)
    slong fmpq_mpoly_factor_get_exp_si(fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_sort(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    int fmpq_mpoly_factor_make_monic(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    int fmpq_mpoly_factor_make_integral(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    int fmpq_mpoly_factor_squarefree(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
    int fmpq_mpoly_factor(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
