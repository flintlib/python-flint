from flint.flintlib.acb cimport acb_t, acb_ptr
from flint.flintlib.acb_poly cimport acb_poly_t
from flint.flintlib.fmpz_poly cimport fmpz_poly_t

cdef extern from "flint/acb_modular.h":
    void acb_modular_theta(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau, long prec)
    void acb_modular_theta_jet(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau, long len, long prec)
    void acb_modular_theta_series(acb_poly_t theta1, acb_poly_t theta2, acb_poly_t theta3, acb_poly_t theta4, const acb_poly_t z, const acb_t tau, long len, long prec)
    void acb_modular_eta(acb_t r, const acb_t tau, long prec)
    void acb_modular_j(acb_t r, const acb_t tau, long prec)
    void acb_modular_lambda(acb_t r, const acb_t tau, long prec)
    void acb_modular_delta(acb_t r, const acb_t tau, long prec)
    void acb_modular_eisenstein(acb_ptr r, const acb_t tau, long len, long prec)
    void acb_modular_elliptic_p(acb_t wp, const acb_t z, const acb_t tau, long prec)
    void acb_modular_elliptic_p_zpx(acb_ptr wp, const acb_t z, const acb_t tau, long len, long prec)
    void acb_modular_elliptic_k(acb_t w, const acb_t m, long prec)
    void acb_modular_elliptic_k_cpx(acb_ptr w, const acb_t m, long len, long prec)
    void acb_modular_elliptic_e(acb_t w, const acb_t m, long prec)
    void acb_modular_hilbert_class_poly(fmpz_poly_t res, long D)
