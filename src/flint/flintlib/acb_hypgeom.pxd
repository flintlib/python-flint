from flint.flintlib.acb cimport acb_t, acb_srcptr, acb_ptr
from flint.flintlib.acb_poly cimport acb_poly_t, acb_poly_struct
from flint.flintlib.mag cimport mag_t

cdef extern from "acb_hypgeom.h":
    void acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, long prec)
    void acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, long prec)
    void acb_hypgeom_bessel_i(acb_t res, const acb_t nu, const acb_t z, long prec)
    void acb_hypgeom_bessel_y(acb_t res, const acb_t nu, const acb_t z, long prec)

    void acb_hypgeom_bessel_k_scaled(acb_t res, const acb_t nu, const acb_t z, long prec)
    void acb_hypgeom_bessel_i_scaled(acb_t res, const acb_t nu, const acb_t z, long prec)

    void acb_hypgeom_erf(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long n, long prec)
    void acb_hypgeom_u_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, long n, long prec)
    void acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t z, long prec)
    void acb_hypgeom_m(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)
    void acb_hypgeom_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)

    long acb_hypgeom_pfq_choose_n(acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long prec)
    void acb_hypgeom_pfq(acb_t res, acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, int regularized, long prec)
    void acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z, int modified, long prec)
    void acb_hypgeom_gamma_upper_asymp(acb_t res, const acb_t s, const acb_t z, int modified, long prec)
    void acb_hypgeom_expint(acb_t res, const acb_t s, const acb_t z, long prec)
    void acb_hypgeom_gamma_lower(acb_t res, const acb_t s, const acb_t z, int modified, long prec)
    void acb_hypgeom_beta_lower(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)
    void acb_hypgeom_erfc(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_erfi(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_pfq_series_direct(acb_poly_t res, const acb_poly_struct * a, long p, const acb_poly_struct * b, long q, const acb_poly_t z, int regularized, long n, long len, long prec)
    void acb_hypgeom_ei(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_si(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_ci(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_shi(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_chi(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_li(acb_t res, const acb_t z, int offset, long prec)
    void acb_hypgeom_2f1(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, long prec)
    void acb_hypgeom_0f1(acb_t res, const acb_t a, const acb_t z, int regularized, long prec)
    void acb_hypgeom_legendre_p(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, long prec)
    void acb_hypgeom_legendre_q(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, long prec)
    void acb_hypgeom_spherical_y(acb_t res, long n, long m, const acb_t theta, const acb_t phi, long prec)
    void acb_hypgeom_jacobi_p(acb_t res, const acb_t n, const acb_t a, const acb_t b, const acb_t z, long prec)
    void acb_hypgeom_gegenbauer_c(acb_t res, const acb_t n, const acb_t m, const acb_t z, long prec)
    void acb_hypgeom_laguerre_l(acb_t res, const acb_t n, const acb_t m, const acb_t z, long prec)
    void acb_hypgeom_hermite_h(acb_t res, const acb_t n, const acb_t z, long prec)
    void acb_hypgeom_chebyshev_t(acb_t res, const acb_t n, const acb_t z, long prec)
    void acb_hypgeom_chebyshev_u(acb_t res, const acb_t n, const acb_t z, long prec)

    void acb_hypgeom_airy_bound(mag_t ai, mag_t aip, mag_t bi, mag_t bip, const acb_t z)
    void acb_hypgeom_airy_asymp(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, long n, long prec)
    void acb_hypgeom_airy_direct(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, long n, long prec)
    void acb_hypgeom_airy(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, long prec)
    void acb_hypgeom_airy_jet(acb_ptr ai, acb_ptr bi, const acb_t z, long len, long prec)
    void _acb_hypgeom_airy_series(acb_ptr ai, acb_ptr ai_prime, acb_ptr bi, acb_ptr bi_prime, acb_srcptr z, long zlen, long len, long prec)
    void acb_hypgeom_airy_series(acb_poly_t ai, acb_poly_t ai_prime, acb_poly_t bi, acb_poly_t bi_prime, const acb_poly_t z, long len, long prec)

    void acb_hypgeom_coulomb(acb_t F, acb_t G, acb_t Hpos, acb_t Hneg, const acb_t l, const acb_t eta, const acb_t z, long prec)
    void acb_hypgeom_coulomb_series(acb_poly_t F, acb_poly_t G, acb_poly_t Hpos, acb_poly_t Hneg, const acb_t l, const acb_t eta, const acb_poly_t z, long len, long prec)

    void acb_hypgeom_erf_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_erfc_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_erfi_series(acb_poly_t res, const acb_poly_t h, long n, long prec)

    void acb_hypgeom_fresnel(acb_t res1, acb_t res2, const acb_t z, int normalized, long prec)
    void acb_hypgeom_fresnel_series(acb_poly_t res1, acb_poly_t res2, const acb_poly_t h, int normalized, long n, long prec)

    void _acb_hypgeom_gamma_upper_series(acb_ptr g, const acb_t s, acb_srcptr h, long hlen, int regularized, long n, long prec)
    void acb_hypgeom_gamma_upper_series(acb_poly_t g, const acb_t s, const acb_poly_t h, int regularized, long n, long prec)

    void _acb_hypgeom_gamma_lower_series(acb_ptr g, const acb_t s, acb_srcptr h, long hlen, int regularized, long n, long prec)
    void acb_hypgeom_gamma_lower_series(acb_poly_t g, const acb_t s, const acb_poly_t h, int regularized, long n, long prec)

    void _acb_hypgeom_beta_lower_series(acb_ptr g, const acb_t s, const acb_t t, acb_srcptr h, long hlen, int regularized, long n, long prec)
    void acb_hypgeom_beta_lower_series(acb_poly_t g, const acb_t s, const acb_t t, const acb_poly_t h, int regularized, long n, long prec)

    void acb_hypgeom_ei_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_si_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_ci_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_shi_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_chi_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_li_series(acb_poly_t res, const acb_poly_t h, int offset, long n, long prec)
