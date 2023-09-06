from flint.flintlib.arb cimport arb_ptr, arb_struct
from flint.flintlib.fmpz_mat cimport fmpz_mat_t
from flint.flintlib.fmpq_mat cimport fmpq_mat_t
from flint.flintlib.mag cimport mag_t
from flint.flintlib.flint cimport ulong
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.arb cimport arb_t
from flint.flintlib.arb_poly cimport arb_poly_t

cdef extern from "arb_mat.h":
    ctypedef struct arb_mat_struct:
        arb_ptr entries
        long r
        long c
        arb_ptr * rows

    ctypedef arb_mat_struct arb_mat_t[1]

    arb_struct * arb_mat_entry(arb_mat_t mat, long i, long j)

    long arb_mat_nrows(const arb_mat_t x)
    long arb_mat_ncols(const arb_mat_t x)

    void arb_mat_init(arb_mat_t mat, long r, long c)
    void arb_mat_clear(arb_mat_t mat)

    void arb_mat_set(arb_mat_t dest, const arb_mat_t src)
    void arb_mat_set_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src)
    void arb_mat_set_fmpq_mat(arb_mat_t dest, const fmpq_mat_t src, long prec)
    void arb_mat_printd(const arb_mat_t mat, long digits)
    int arb_mat_equal(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_overlaps(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_contains(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_contains_fmpq_mat(const arb_mat_t mat1, const fmpq_mat_t mat2)
    int arb_mat_contains_fmpz_mat(const arb_mat_t mat1, const fmpz_mat_t mat2)

    void arb_mat_zero(arb_mat_t mat)
    void arb_mat_one(arb_mat_t mat)

    void arb_mat_bound_inf_norm(mag_t b, const arb_mat_t A)

    void arb_mat_neg(arb_mat_t dest, const arb_mat_t src)
    void arb_mat_add(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    void arb_mat_sub(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    void arb_mat_mul(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    void arb_mat_pow_ui(arb_mat_t B, const arb_mat_t A, ulong exp, long prec)

    void arb_mat_scalar_mul_2exp_si(arb_mat_t B, const arb_mat_t A, long c)
    void arb_mat_scalar_addmul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)
    void arb_mat_scalar_mul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)
    void arb_mat_scalar_div_si(arb_mat_t B, const arb_mat_t A, long c, long prec)
    void arb_mat_scalar_addmul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)
    void arb_mat_scalar_mul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)
    void arb_mat_scalar_div_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)
    void arb_mat_scalar_addmul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
    void arb_mat_scalar_mul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
    void arb_mat_scalar_div_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)

    int arb_mat_lu(long * P, arb_mat_t LU, const arb_mat_t A, long prec)
    void arb_mat_solve_lu_precomp(arb_mat_t X, const long * perm, const arb_mat_t A, const arb_mat_t B, long prec)
    int arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
    int arb_mat_solve_lu(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
    int arb_mat_solve_precond(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
    int arb_mat_inv(arb_mat_t X, const arb_mat_t A, long prec)
    void arb_mat_det(arb_t det, const arb_mat_t A, long prec)

    void arb_mat_exp(arb_mat_t B, const arb_mat_t A, long prec)

    void _arb_mat_charpoly(arb_ptr cp, const arb_mat_t mat, long prec)
    void arb_mat_charpoly(arb_poly_t cp, const arb_mat_t mat, long prec)

    void arb_mat_transpose(arb_mat_t B, const arb_mat_t A)

    void arb_mat_trace(arb_t trace, const arb_mat_t mat, long prec)
    void arb_mat_ones(arb_mat_t mat)
    void arb_mat_hilbert(arb_mat_t mat, long prec)
    void arb_mat_pascal(arb_mat_t mat, int triangular, long prec)
    void arb_mat_stirling(arb_mat_t mat, int kind, long prec)
    void arb_mat_dct(arb_mat_t mat, int type, long prec)

    void arb_mat_get_mid(arb_mat_t B, const arb_mat_t A)

    int arb_mat_eq(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_ne(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_equal(const arb_mat_t mat1, const arb_mat_t mat2)

    void arb_mat_frobenius_norm(arb_t res, const arb_mat_t A, long prec)

    int arb_mat_approx_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
