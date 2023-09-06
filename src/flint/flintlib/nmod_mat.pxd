from flint.flintlib.flint cimport mp_limb_t, flint_rand_t, mp_ptr
from flint.flintlib.flint cimport mp_srcptr
from flint.flintlib.nmod_vec cimport nmod_t
from flint.flintlib.nmod_poly cimport nmod_poly_t

cdef extern from "flint/nmod_mat.h":
    ctypedef struct nmod_mat_struct:
        mp_limb_t * entries
        long r
        long c
        mp_limb_t ** rows
        nmod_t mod
    ctypedef nmod_mat_struct nmod_mat_t[1]
    mp_limb_t nmod_mat_entry(nmod_mat_t mat, long i, long j)
    long nmod_mat_nrows(nmod_mat_t mat)
    long nmod_mat_ncols(nmod_mat_t mat)
    void _nmod_mat_set_mod(nmod_mat_t mat, mp_limb_t n)
    void nmod_mat_init(nmod_mat_t mat, long rows, long cols, mp_limb_t n)
    void nmod_mat_init_set(nmod_mat_t mat, nmod_mat_t src)
    void nmod_mat_clear(nmod_mat_t mat)
    void nmod_mat_randtest(nmod_mat_t mat, flint_rand_t state)
    void nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state)
    void nmod_mat_randrank(nmod_mat_t, flint_rand_t state, long rank)
    void nmod_mat_randops(nmod_mat_t mat, long count, flint_rand_t state)
    void nmod_mat_randtril(nmod_mat_t mat, flint_rand_t state, int unit)
    void nmod_mat_randtriu(nmod_mat_t mat, flint_rand_t state, int unit)
    void nmod_mat_print_pretty(nmod_mat_t mat)
    int nmod_mat_equal(nmod_mat_t mat1, nmod_mat_t mat2)
    int nmod_mat_is_zero(nmod_mat_t mat)
    int nmod_mat_is_empty(nmod_mat_t mat)
    int nmod_mat_is_square(nmod_mat_t mat)
    void nmod_mat_zero(nmod_mat_t mat)
    void nmod_mat_set(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_transpose(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_add(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_sub(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_neg(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_scalar_mul(nmod_mat_t B, nmod_mat_t A, mp_limb_t c)
    void nmod_mat_mul(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_mul_classical(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_mul_strassen(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_addmul(nmod_mat_t D, nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_submul(nmod_mat_t D, nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    mp_limb_t nmod_mat_det(nmod_mat_t A)
    long nmod_mat_rank(nmod_mat_t A)
    int nmod_mat_inv(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_solve_tril(nmod_mat_t X, nmod_mat_t L, nmod_mat_t B, int unit)
    void nmod_mat_solve_triu(nmod_mat_t X, nmod_mat_t U, nmod_mat_t B, int unit)
    long nmod_mat_lu(long * P, nmod_mat_t A, int rank_check)
    int nmod_mat_solve(nmod_mat_t X, nmod_mat_t A, nmod_mat_t B)
    int nmod_mat_solve_vec(mp_ptr x, nmod_mat_t A, mp_srcptr b)
    long nmod_mat_rref(nmod_mat_t A)
    long nmod_mat_nullspace(nmod_mat_t X, nmod_mat_t A)
