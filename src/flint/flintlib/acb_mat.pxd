from flint.flintlib.flint cimport ulong, flint_rand_t
from flint.flintlib.fmpz_mat cimport fmpz_mat_t
from flint.flintlib.fmpq_mat cimport fmpq_mat_t
from flint.flintlib.mag cimport mag_t
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.acb_poly cimport acb_poly_t
from flint.flintlib.arb cimport arb_t
from flint.flintlib.acb cimport  acb_ptr, acb_struct, acb_t, acb_srcptr

cdef extern from "acb_mat.h":
    ctypedef struct acb_mat_struct:
        acb_ptr entries
        long r
        long c
        acb_ptr * rows

    ctypedef acb_mat_struct acb_mat_t[1]

    acb_struct * acb_mat_entry(acb_mat_t mat, long i, long j)

    long acb_mat_nrows(const acb_mat_t x)
    long acb_mat_ncols(const acb_mat_t x)

    void acb_mat_init(acb_mat_t mat, long r, long c)
    void acb_mat_clear(acb_mat_t mat)

    void acb_mat_set(acb_mat_t dest, const acb_mat_t src)
    void acb_mat_set_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src)
    void acb_mat_set_fmpq_mat(acb_mat_t dest, const fmpq_mat_t src, long prec)
    void acb_mat_printd(const acb_mat_t mat, long digits)
    int acb_mat_equal(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_overlaps(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_contains(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_contains_fmpq_mat(const acb_mat_t mat1, const fmpq_mat_t mat2)
    int acb_mat_contains_fmpz_mat(const acb_mat_t mat1, const fmpz_mat_t mat2)

    void acb_mat_zero(acb_mat_t mat)
    void acb_mat_one(acb_mat_t mat)

    void acb_mat_bound_inf_norm(mag_t b, const acb_mat_t A)

    void acb_mat_neg(acb_mat_t dest, const acb_mat_t src)
    void acb_mat_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_pow_ui(acb_mat_t B, const acb_mat_t A, ulong exp, long prec)

    void acb_mat_scalar_mul_2exp_si(acb_mat_t B, const acb_mat_t A, long c)
    void acb_mat_scalar_addmul_si(acb_mat_t B, const acb_mat_t A, long c, long prec)
    void acb_mat_scalar_mul_si(acb_mat_t B, const acb_mat_t A, long c, long prec)
    void acb_mat_scalar_div_si(acb_mat_t B, const acb_mat_t A, long c, long prec)
    void acb_mat_scalar_addmul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)
    void acb_mat_scalar_mul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)
    void acb_mat_scalar_div_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)
    void acb_mat_scalar_addmul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)
    void acb_mat_scalar_mul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)
    void acb_mat_scalar_div_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)

    int acb_mat_lu(long * P, acb_mat_t LU, const acb_mat_t A, long prec)
    void acb_mat_solve_lu_precomp(acb_mat_t X, const long * perm, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_solve_lu(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_solve_precond(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_inv(acb_mat_t X, const acb_mat_t A, long prec)
    void acb_mat_det(acb_t det, const acb_mat_t A, long prec)

    void acb_mat_exp(acb_mat_t B, const acb_mat_t A, long prec)

    void _acb_mat_charpoly(acb_ptr cp, const acb_mat_t mat, long prec)
    void acb_mat_charpoly(acb_poly_t cp, const acb_mat_t mat, long prec)


    void acb_mat_conjugate(acb_mat_t mat1, const acb_mat_t mat2)
    void acb_mat_transpose(acb_mat_t B, const acb_mat_t A)
    void acb_mat_trace(acb_t trace, const acb_mat_t mat, long prec)
    void acb_mat_get_mid(acb_mat_t B, const acb_mat_t A)

    void acb_mat_dft(acb_mat_t res, int kind, long prec)

    void acb_mat_frobenius_norm(arb_t res, const acb_mat_t A, long prec)

    int acb_mat_eq(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_ne(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_equal(const acb_mat_t mat1, const acb_mat_t mat2)

    int acb_mat_approx_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)

    void acb_mat_randtest_eig(acb_mat_t A, flint_rand_t state, acb_srcptr E, long prec)
    int acb_mat_approx_eig_qr(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, const mag_t tol, long maxiter, long prec)
    void acb_mat_eig_global_enclosure(mag_t eps, const acb_mat_t A, acb_srcptr E, const acb_mat_t R, long prec)

    int acb_mat_eig_simple_rump(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)
    int acb_mat_eig_simple_vdhoeven_mourrain(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)
    int acb_mat_eig_simple(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)

    int acb_mat_eig_multiple_rump(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)
    int acb_mat_eig_multiple(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)
