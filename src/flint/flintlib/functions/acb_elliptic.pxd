from flint.flintlib.types.acb cimport acb_poly_t, acb_ptr, acb_srcptr, acb_t
from flint.flintlib.types.flint cimport slong



cdef extern from "flint/acb_elliptic.h":
    void acb_elliptic_k(acb_t res, const acb_t m, slong prec)
    void acb_elliptic_k_jet(acb_ptr res, const acb_t m, slong len, slong prec)
    void _acb_elliptic_k_series(acb_ptr res, acb_srcptr m, slong mlen, slong len, slong prec)
    void acb_elliptic_k_series(acb_poly_t res, const acb_poly_t m, slong len, slong prec)
    void acb_elliptic_e(acb_t res, const acb_t m, slong prec)
    void acb_elliptic_pi(acb_t res, const acb_t n, const acb_t m, slong prec)
    void acb_elliptic_f(acb_t res, const acb_t phi, const acb_t m, int pi, slong prec)
    void acb_elliptic_e_inc(acb_t res, const acb_t phi, const acb_t m, int pi, slong prec)
    void acb_elliptic_pi_inc(acb_t res, const acb_t n, const acb_t phi, const acb_t m, int pi, slong prec)
    void acb_elliptic_rf(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, slong prec)
    void acb_elliptic_rg(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, slong prec)
    void acb_elliptic_rj(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, slong prec)
    void acb_elliptic_rj_carlson(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, slong prec)
    void acb_elliptic_rj_integration(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, slong prec)
    void acb_elliptic_rc1(acb_t res, const acb_t x, slong prec)
    void acb_elliptic_p(acb_t res, const acb_t z, const acb_t tau, slong prec)
    void acb_elliptic_p_prime(acb_t res, const acb_t z, const acb_t tau, slong prec)
    void acb_elliptic_p_jet(acb_ptr res, const acb_t z, const acb_t tau, slong len, slong prec)
    void _acb_elliptic_p_series(acb_ptr res, acb_srcptr z, slong zlen, const acb_t tau, slong len, slong prec)
    void acb_elliptic_p_series(acb_poly_t res, const acb_poly_t z, const acb_t tau, slong len, slong prec)
    void acb_elliptic_invariants(acb_t g2, acb_t g3, const acb_t tau, slong prec)
    void acb_elliptic_roots(acb_t e1, acb_t e2, acb_t e3, const acb_t tau, slong prec)
    void acb_elliptic_inv_p(acb_t res, const acb_t z, const acb_t tau, slong prec)
    void acb_elliptic_zeta(acb_t res, const acb_t z, const acb_t tau, slong prec)
    void acb_elliptic_sigma(acb_t res, const acb_t z, const acb_t tau, slong prec)
