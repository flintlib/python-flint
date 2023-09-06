from flint.flintlib.arb cimport arb_t, arb_srcptr
from flint.flintlib.arb_poly cimport arb_poly_t
from flint.flintlib.fmpz cimport fmpz_t
from flint._flint cimport ulong

cdef extern from "flint/arb_hypgeom.h":
    void arb_hypgeom_pfq(arb_t res, arb_srcptr a, long p, arb_srcptr b, long q, const arb_t z, int regularized, long prec)
    void arb_hypgeom_0f1(arb_t res, const arb_t a, const arb_t z, int regularized, long prec)
    void arb_hypgeom_m(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)
    void arb_hypgeom_1f1(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)
    void arb_hypgeom_u(arb_t res, const arb_t a, const arb_t b, const arb_t z, long prec)
    void arb_hypgeom_2f1(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, long prec)

    void arb_hypgeom_erf(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_erf_series(arb_poly_t g, const arb_poly_t h, long len, long prec)
    void arb_hypgeom_erfc(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_erfc_series(arb_poly_t g, const arb_poly_t h, long len, long prec)
    void arb_hypgeom_erfi(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_erfi_series(arb_poly_t g, const arb_poly_t h, long len, long prec)
    void arb_hypgeom_fresnel(arb_t res1, arb_t res2, const arb_t z, int normalized, long prec)
    void arb_hypgeom_fresnel_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, int normalized, long len, long prec)

    void arb_hypgeom_ei(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_si(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_ci(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_shi(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_chi(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_li(arb_t res, const arb_t z, int offset, long prec)
    void arb_hypgeom_ei_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_si_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_ci_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_shi_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_chi_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_li_series(arb_poly_t res, const arb_poly_t h, int offset, long n, long prec)

    void arb_hypgeom_bessel_j(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_bessel_k(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_bessel_i(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_bessel_y(arb_t res, const arb_t nu, const arb_t z, long prec)

    void arb_hypgeom_bessel_k_scaled(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_bessel_i_scaled(arb_t res, const arb_t nu, const arb_t z, long prec)

    void arb_hypgeom_airy(arb_t ai, arb_t aip, arb_t bi, arb_t bip, const arb_t z, long prec)
    void arb_hypgeom_airy_series(arb_poly_t ai, arb_poly_t ai_prime, arb_poly_t bi, arb_poly_t bi_prime, const arb_poly_t z, long len, long prec)
    void arb_hypgeom_airy_zero(arb_t ai, arb_t aip, arb_t bi, arb_t bip, const fmpz_t n, long prec)

    void arb_hypgeom_coulomb(arb_t F, arb_t G, const arb_t l, const arb_t eta, const arb_t z, long prec)
    void arb_hypgeom_coulomb_series(arb_poly_t F, arb_poly_t G, const arb_t l, const arb_t eta, const arb_poly_t z, long len, long prec)

    void arb_hypgeom_expint(arb_t res, const arb_t s, const arb_t z, long prec)
    void arb_hypgeom_gamma_upper(arb_t res, const arb_t s, const arb_t z, int modified, long prec)
    void arb_hypgeom_gamma_lower(arb_t res, const arb_t s, const arb_t z, int modified, long prec)
    void arb_hypgeom_beta_lower(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)

    void arb_hypgeom_gamma_upper_series(arb_poly_t g, const arb_t s, const arb_poly_t h, int regularized, long n, long prec)
    void arb_hypgeom_gamma_lower_series(arb_poly_t g, const arb_t s, const arb_poly_t h, int regularized, long n, long prec)
    void arb_hypgeom_beta_lower_series(arb_poly_t g, const arb_t s, const arb_t t, const arb_poly_t h, int regularized, long n, long prec)

    void arb_hypgeom_chebyshev_t(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_chebyshev_u(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_jacobi_p(arb_t res, const arb_t n, const arb_t a, const arb_t b, const arb_t z, long prec)
    void arb_hypgeom_gegenbauer_c(arb_t res, const arb_t n, const arb_t m, const arb_t z, long prec)
    void arb_hypgeom_laguerre_l(arb_t res, const arb_t n, const arb_t m, const arb_t z, long prec)
    void arb_hypgeom_hermite_h(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_legendre_p(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, long prec)
    void arb_hypgeom_legendre_q(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, long prec)

    void arb_hypgeom_legendre_p_ui_root(arb_t res, arb_t weight, ulong n, ulong k, long prec)
