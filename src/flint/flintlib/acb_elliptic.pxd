from flint.flintlib.acb cimport acb_t

cdef extern from "acb_elliptic.h":
    void acb_elliptic_rf(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, long prec)
    void acb_elliptic_rj(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, long prec)
    void acb_elliptic_rg(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, long prec)
    void acb_elliptic_f(acb_t res, const acb_t phi, const acb_t m, int times_pi, long prec)
    void acb_elliptic_e_inc(acb_t res, const acb_t phi, const acb_t m, int times_pi, long prec)
    void acb_elliptic_pi(acb_t res, const acb_t n, const acb_t m, long prec)
    void acb_elliptic_pi_inc(acb_t res, const acb_t n, const acb_t phi, const acb_t m, int times_pi, long prec)
    void acb_elliptic_p(acb_t res, const acb_t z, const acb_t tau, long prec)
    void acb_elliptic_zeta(acb_t res, const acb_t z, const acb_t tau, long prec)
    void acb_elliptic_sigma(acb_t res, const acb_t z, const acb_t tau, long prec)
    void acb_elliptic_roots(acb_t e1, acb_t e2, acb_t e3, const acb_t tau, long prec)
    void acb_elliptic_invariants(acb_t g2, acb_t g3, const acb_t tau, long prec)
    void acb_elliptic_inv_p(acb_t res, const acb_t z, const acb_t tau, long prec)
