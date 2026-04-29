from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.fmpz_mod cimport fmpz_mod_mpoly_ctx_t, fmpz_mod_mpoly_t
from flint.flintlib.types.gr cimport gr_ctx_t


cdef extern from *:
    """
    #include "flint/gr.h"
    #if __FLINT_RELEASE >= 30400 /* Flint 3.4.0 or later */
    #include "flint/gr_series.h"
    #endif

    #if __FLINT_RELEASE >= 30400 /* Flint 3.4.0 or later */

    #define compat_gr_ctx_init_gr_series(...) gr_series_ctx_init(__VA_ARGS__)
    #define compat_gr_ctx_init_series_mod_gr_poly(...) gr_series_mod_ctx_init(__VA_ARGS__)

    #else

    #define compat_gr_ctx_init_gr_series(...) gr_ctx_init_gr_series(__VA_ARGS__)
    #define compat_gr_ctx_init_series_mod_gr_poly(...) gr_ctx_init_series_mod_gr_poly(__VA_ARGS__)

    #endif

    #if __FLINT_RELEASE < 30200 /* Flint < 3.2.0 */

    #define compat_fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(...) (void)0

    #else

    #define compat_fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(...) fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(__VA_ARGS__)

    #endif
    """
    void compat_gr_ctx_init_gr_series(gr_ctx_t ctx, gr_ctx_t base_ring, slong prec)
    void compat_gr_ctx_init_series_mod_gr_poly(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)
    void compat_fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const slong * c, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)
