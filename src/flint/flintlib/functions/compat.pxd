from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.fmpz_mod cimport fmpz_mod_mpoly_ctx_t, fmpz_mod_mpoly_t


cdef extern from *:
    """
    #if __FLINT_RELEASE < 30200 /* Flint < 3.2.0 */

    #define compat_fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(...) (void)0

    #else

    #define compat_fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(...) fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(__VA_ARGS__)

    #endif
    """
    void compat_fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const slong * c, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)
