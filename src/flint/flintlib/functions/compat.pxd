from flint.flintlib.types.gr cimport gr_ctx_t, gr_ptr, gr_srcptr, gr_vec_t


cdef extern from *:
    """
    #include "flint/gr.h"
    #if __FLINT_RELEASE >= 30600 /* Flint 3.6.0 or later */

    #define compat_gr_factor(c, factors, exponents, x, flags, ctx) gr_factor(c, factors, (fmpz_vec_struct *) exponents, x, flags, ctx)

    #else

    #define compat_gr_factor(c, factors, exponents, x, flags, ctx) gr_factor(c, factors, exponents, x, flags, ctx)

    #endif
    """
    int compat_gr_factor(gr_ptr c, gr_vec_t factors, gr_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx)
