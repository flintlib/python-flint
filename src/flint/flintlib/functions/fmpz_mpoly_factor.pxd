from flint.flintlib.types.flint cimport fmpz_t, slong
from flint.flintlib.types.fmpq cimport fmpq_t
from flint.flintlib.types.fmpz cimport fmpz_mpoly_ctx_t, fmpz_mpoly_factor_t, fmpz_mpoly_t



cdef extern from "flint/fmpz_mpoly_factor.h":
    void fmpz_mpoly_factor_init(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_swap(fmpz_mpoly_factor_t f, fmpz_mpoly_factor_t g, const fmpz_mpoly_ctx_t ctx)
    slong fmpz_mpoly_factor_length(const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_get_constant_fmpz(fmpz_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_get_base(fmpz_mpoly_t B, const fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_swap_base(fmpz_mpoly_t B, fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t ctx)
    slong fmpz_mpoly_factor_get_exp_si(fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_factor_sort(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_factor_squarefree(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_factor(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
