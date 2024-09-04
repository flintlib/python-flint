from flint.flintlib.types.acb cimport acb_ptr, acb_srcptr
from flint.flintlib.types.flint cimport slong

# unknown type acb_dft_bluestein_t
# unknown type acb_dft_crt_t
# unknown type acb_dft_cyc_t
# unknown type acb_dft_naive_t
# unknown type acb_dft_pre_t
# unknown type acb_dft_prod_t
# unknown type acb_dft_rad2_t


cdef extern from "flint/acb_dft.h":
    void acb_dft(acb_ptr w, acb_srcptr v, slong n, slong prec)
    void acb_dft_inverse(acb_ptr w, acb_srcptr v, slong n, slong prec)
    # void acb_dft_precomp_init(acb_dft_pre_t pre, slong len, slong prec)
    # void acb_dft_precomp_clear(acb_dft_pre_t pre)
    # void acb_dft_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)
    # void acb_dft_inverse_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)
    void acb_dirichlet_dft_prod(acb_ptr w, acb_srcptr v, slong * cyc, slong num, slong prec)
    # void acb_dft_prod_init(acb_dft_prod_t t, slong * cyc, slong num, slong prec)
    # void acb_dft_prod_clear(acb_dft_prod_t t)
    # void acb_dirichlet_dft_prod_precomp(acb_ptr w, acb_srcptr v, const acb_dft_prod_t prod, slong prec)
    void acb_dft_convol_naive(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
    void acb_dft_convol_rad2(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
    void acb_dft_convol(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
    void acb_dft_naive(acb_ptr w, acb_srcptr v, slong n, slong prec)
    # void acb_dft_naive_init(acb_dft_naive_t t, slong len, slong prec)
    # void acb_dft_naive_clear(acb_dft_naive_t t)
    # void acb_dft_naive_precomp(acb_ptr w, acb_srcptr v, const acb_dft_naive_t t, slong prec)
    void acb_dft_crt(acb_ptr w, acb_srcptr v, slong n, slong prec)
    # void acb_dft_crt_init(acb_dft_crt_t t, slong len, slong prec)
    # void acb_dft_crt_clear(acb_dft_crt_t t)
    # void acb_dft_crt_precomp(acb_ptr w, acb_srcptr v, const acb_dft_crt_t t, slong prec)
    void acb_dft_cyc(acb_ptr w, acb_srcptr v, slong n, slong prec)
    # void acb_dft_cyc_init(acb_dft_cyc_t t, slong len, slong prec)
    # void acb_dft_cyc_clear(acb_dft_cyc_t t)
    # void acb_dft_cyc_precomp(acb_ptr w, acb_srcptr v, const acb_dft_cyc_t t, slong prec)
    void acb_dft_rad2(acb_ptr w, acb_srcptr v, int e, slong prec)
    void acb_dft_inverse_rad2(acb_ptr w, acb_srcptr v, int e, slong prec)
    # void acb_dft_rad2_init(acb_dft_rad2_t t, int e, slong prec)
    # void acb_dft_rad2_clear(acb_dft_rad2_t t)
    # void acb_dft_rad2_precomp(acb_ptr w, acb_srcptr v, const acb_dft_rad2_t t, slong prec)
    void acb_dft_bluestein(acb_ptr w, acb_srcptr v, slong n, slong prec)
    # void acb_dft_bluestein_init(acb_dft_bluestein_t t, slong len, slong prec)
    # void acb_dft_bluestein_clear(acb_dft_bluestein_t t)
    # void acb_dft_bluestein_precomp(acb_ptr w, acb_srcptr v, const acb_dft_bluestein_t t, slong prec)
