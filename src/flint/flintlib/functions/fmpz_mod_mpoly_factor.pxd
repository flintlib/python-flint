from flint.flintlib.types.flint cimport fmpz_t, slong
from flint.flintlib.types.fmpz_mod cimport fmpz_mod_mpoly_ctx_t, fmpz_mod_mpoly_factor_t, fmpz_mod_mpoly_t



cdef extern from "flint/fmpz_mod_mpoly_factor.h":
    void fmpz_mod_mpoly_factor_init(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_factor_clear(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_factor_swap(fmpz_mod_mpoly_factor_t f, fmpz_mod_mpoly_factor_t g, const fmpz_mod_mpoly_ctx_t ctx)
    slong fmpz_mod_mpoly_factor_length(const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_factor_get_constant_fmpz(fmpz_t c, const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_factor_get_base(fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_factor_swap_base(fmpz_mod_mpoly_t B, fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    slong fmpz_mod_mpoly_factor_get_exp_si(fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)
    void fmpz_mod_mpoly_factor_sort(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
    int fmpz_mod_mpoly_factor_squarefree(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
    int fmpz_mod_mpoly_factor(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
