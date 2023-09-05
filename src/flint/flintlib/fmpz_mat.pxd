from flint._flint cimport flint_rand_t, mp_bitcnt_t, ulong
from flint.flintlib.fmpz cimport fmpz_struct, fmpz_t
from flint.flintlib.fmpz_poly cimport fmpz_poly_t
from flint.flintlib.nmod_mat cimport nmod_mat_t

cdef extern from "flint/fmpz_mat.h":
    ctypedef struct fmpz_mat_struct:
        fmpz_struct * entries
        long r
        long c
        fmpz_struct ** rows
    ctypedef fmpz_mat_struct fmpz_mat_t[1]
    fmpz_struct * fmpz_mat_entry(fmpz_mat_t mat, long i, long j)
    long fmpz_mat_nrows(fmpz_mat_t mat)
    long fmpz_mat_ncols(fmpz_mat_t mat)
    void fmpz_mat_init(fmpz_mat_t mat, long rows, long cols)
    void fmpz_mat_init_set(fmpz_mat_t mat,  fmpz_mat_t src)
    void fmpz_mat_swap(fmpz_mat_t mat1, fmpz_mat_t mat2)
    void fmpz_mat_set(fmpz_mat_t mat1,  fmpz_mat_t mat2)
    void fmpz_mat_clear(fmpz_mat_t mat)
    int fmpz_mat_equal(fmpz_mat_t mat1, fmpz_mat_t mat2)
    int fmpz_mat_is_zero( fmpz_mat_t mat)
    int fmpz_mat_is_empty( fmpz_mat_t mat)
    int fmpz_mat_is_square( fmpz_mat_t mat)
    void fmpz_mat_zero(fmpz_mat_t mat)
    void fmpz_mat_one(fmpz_mat_t mat)
    void fmpz_mat_get_nmod_mat(nmod_mat_t mat, fmpz_mat_t mat2)
    void fmpz_mat_randbits(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpz_mat_randtest(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpz_mat_randtest_unsigned(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpz_mat_randrank(fmpz_mat_t mat, flint_rand_t state, long rank, mp_bitcnt_t bits)
    int fmpz_mat_print_pretty( fmpz_mat_t mat)
    void fmpz_mat_transpose(fmpz_mat_t B,  fmpz_mat_t A)
    void fmpz_mat_add(fmpz_mat_t C,  fmpz_mat_t A,  fmpz_mat_t B)
    void fmpz_mat_sub(fmpz_mat_t C,  fmpz_mat_t A,  fmpz_mat_t B)
    void fmpz_mat_neg(fmpz_mat_t B,  fmpz_mat_t A)
    void fmpz_mat_scalar_mul_fmpz(fmpz_mat_t B,  fmpz_mat_t A,  fmpz_t c)
    void fmpz_mat_scalar_mul_si(fmpz_mat_t B,  fmpz_mat_t A, long c)
    void fmpz_mat_scalar_mul_ui(fmpz_mat_t B,  fmpz_mat_t A, ulong c)
    void fmpz_mat_mul(fmpz_mat_t C,  fmpz_mat_t A,  fmpz_mat_t B)
    void fmpz_mat_det(fmpz_t det,  fmpz_mat_t A)
    long fmpz_mat_rank(fmpz_mat_t A)
    long fmpz_mat_rref(fmpz_mat_t R, fmpz_t den, fmpz_mat_t A)
    void fmpz_mat_inv(fmpz_mat_t B, fmpz_t den, fmpz_mat_t A)
    int fmpz_mat_solve(fmpz_mat_t X, fmpz_t den, fmpz_mat_t A, fmpz_mat_t B)
    long fmpz_mat_nullspace(fmpz_mat_t res, fmpz_mat_t mat)
    void fmpz_mat_pow(fmpz_mat_t A, fmpz_mat_t B, ulong e)
    int fmpz_mat_is_hadamard(const fmpz_mat_t A)
    int fmpz_mat_hadamard(fmpz_mat_t A)

    void fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A)
    void fmpz_mat_hnf_transform(fmpz_mat_t H, fmpz_mat_t U, const  fmpz_mat_t A)
    int fmpz_mat_is_in_hnf(const fmpz_mat_t A)
    void fmpz_mat_snf(fmpz_mat_t S, const fmpz_mat_t A)
    int fmpz_mat_is_in_snf(const fmpz_mat_t A)

    void fmpz_mat_charpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
    void fmpz_mat_minpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
