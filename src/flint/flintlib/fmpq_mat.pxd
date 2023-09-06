from flint.flintlib.flint cimport flint_rand_t, mp_bitcnt_t
from flint.flintlib.fmpz cimport fmpz_struct, fmpz_t
from flint.flintlib.fmpq cimport fmpq_struct, fmpq_t
from flint.flintlib.fmpz_mat cimport fmpz_mat_t
from flint.flintlib.fmpq_poly cimport fmpq_poly_t

cdef extern from "flint/fmpq_mat.h":
    ctypedef struct fmpq_mat_struct:
        fmpq_struct * entries
        long r
        long c
        fmpq_struct ** rows
    ctypedef fmpq_mat_struct fmpq_mat_t[1]
    fmpq_struct * fmpq_mat_entry(fmpq_mat_t mat, long i, long j)
    fmpz_struct * fmpq_mat_entry_num(fmpq_mat_t mat, long i, long j)
    fmpz_struct * fmpq_mat_entry_den(fmpq_mat_t mat, long i, long j)
    long fmpq_mat_nrows(fmpq_mat_t mat)
    long fmpq_mat_ncols(fmpq_mat_t mat)
    void fmpq_mat_init(fmpq_mat_t mat, long rows, long cols)
    void fmpq_mat_clear(fmpq_mat_t mat)
    void fmpq_mat_print(fmpq_mat_t mat)
    void fmpq_mat_randbits(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_mat_randtest(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_mat_hilbert_matrix(fmpq_mat_t mat)
    void fmpq_mat_set(fmpq_mat_t dest, fmpq_mat_t src)
    void fmpq_mat_zero(fmpq_mat_t mat)
    void fmpq_mat_one(fmpq_mat_t mat)
    void fmpq_mat_add(fmpq_mat_t mat, fmpq_mat_t mat1, fmpq_mat_t mat2)
    void fmpq_mat_sub(fmpq_mat_t mat, fmpq_mat_t mat1, fmpq_mat_t mat2)
    void fmpq_mat_neg(fmpq_mat_t rop, fmpq_mat_t op)
    void fmpq_mat_scalar_mul_fmpz(fmpq_mat_t rop, fmpq_mat_t op, fmpz_t x)
    void fmpq_mat_scalar_div_fmpz(fmpq_mat_t rop, fmpq_mat_t op, fmpz_t x)
    int fmpq_mat_equal(fmpq_mat_t mat1, fmpq_mat_t mat2)
    int fmpq_mat_is_integral(fmpq_mat_t mat)
    int fmpq_mat_is_zero(fmpq_mat_t mat)
    int fmpq_mat_is_empty(fmpq_mat_t mat)
    int fmpq_mat_is_square(fmpq_mat_t mat)
    int fmpq_mat_get_fmpz_mat(fmpz_mat_t dest, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_entrywise(fmpz_mat_t num, fmpz_mat_t den, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t num, fmpz_t den, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_t num, fmpz_struct * den, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_colwise(fmpz_mat_t num, fmpz_struct * den, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_rowwise_2(fmpz_mat_t num, fmpz_mat_t num2, fmpz_struct * den, fmpq_mat_t mat, fmpq_mat_t mat2)
    void fmpq_mat_get_fmpz_mat_mod_fmpz(fmpz_mat_t dest, fmpq_mat_t mat, fmpz_t mod)
    void fmpq_mat_set_fmpz_mat(fmpq_mat_t dest, fmpz_mat_t src)
    void fmpq_mat_set_fmpz_mat_div_fmpz(fmpq_mat_t X, fmpz_mat_t Xmod, fmpz_t div)
    int fmpq_mat_set_fmpz_mat_mod_fmpz(fmpq_mat_t X, fmpz_mat_t Xmod, fmpz_t mod)
    void fmpq_mat_mul(fmpq_mat_t C, fmpq_mat_t A, fmpq_mat_t B)
    void fmpq_mat_mul_fmpz_mat(fmpq_mat_t C, fmpq_mat_t A, fmpz_mat_t B)
    void fmpq_mat_mul_r_fmpz_mat(fmpq_mat_t C, fmpz_mat_t A, fmpq_mat_t B)
    void fmpq_mat_det(fmpq_t det, fmpq_mat_t mat)
    int fmpq_mat_solve_fraction_free(fmpq_mat_t X, fmpq_mat_t A, fmpq_mat_t B)
    int fmpq_mat_solve_dixon(fmpq_mat_t X, fmpq_mat_t A, fmpq_mat_t B)
    int fmpq_mat_solve_fmpz_mat(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
    int fmpq_mat_inv(fmpq_mat_t B, fmpq_mat_t A)
    long fmpq_mat_rref(fmpq_mat_t B, fmpq_mat_t A)
    void fmpq_mat_transpose(fmpq_mat_t B, fmpq_mat_t A)

    void fmpq_mat_charpoly(fmpq_poly_t cp, const fmpq_mat_t mat)
    void fmpq_mat_minpoly(fmpq_poly_t cp, const fmpq_mat_t mat)
