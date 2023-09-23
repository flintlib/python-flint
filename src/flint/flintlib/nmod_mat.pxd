from flint.flintlib.flint cimport mp_limb_t, flint_rand_t, mp_ptr
from flint.flintlib.flint cimport mp_srcptr, slong, ulong
from flint.flintlib.nmod cimport nmod_t
from flint.flintlib.nmod_poly cimport nmod_poly_t
from flint.flintlib.fmpz cimport fmpz_t

cdef extern from "flint/nmod_mat.h":
    ctypedef struct nmod_mat_struct:
        mp_limb_t * entries
        long r
        long c
        mp_limb_t ** rows
        nmod_t mod
    ctypedef nmod_mat_struct nmod_mat_t[1]

# from here on is parsed
    void nmod_mat_init(nmod_mat_t mat, slong rows, slong cols, mp_limb_t n)
    void nmod_mat_init_set(nmod_mat_t mat, nmod_mat_t src)
    void nmod_mat_clear(nmod_mat_t mat)
    void nmod_mat_set(nmod_mat_t mat, nmod_mat_t src)
    void nmod_mat_swap(nmod_mat_t mat1, nmod_mat_t mat2)
    void nmod_mat_swap_entrywise(nmod_mat_t mat1, nmod_mat_t mat2)
    mp_limb_t nmod_mat_entry(nmod_mat_t mat, long i, long j)
    # MACRO nmod_mat_entry(nmod_mat_t mat, slong i, slong j)
    mp_limb_t nmod_mat_get_entry(const nmod_mat_t mat, slong i, slong j)
    mp_limb_t * nmod_mat_entry_ptr(const nmod_mat_t mat, slong i, slong j)
    void nmod_mat_set_entry(nmod_mat_t mat, slong i, slong j, mp_limb_t x)
    slong nmod_mat_nrows(nmod_mat_t mat)
    slong nmod_mat_ncols(nmod_mat_t mat)
    void nmod_mat_zero(nmod_mat_t mat)
    int nmod_mat_is_zero(nmod_mat_t mat)
    void nmod_mat_window_init(nmod_mat_t window, const nmod_mat_t mat, slong r1, slong c1, slong r2, slong c2)
    void nmod_mat_window_clear(nmod_mat_t window)
    void nmod_mat_concat_vertical(nmod_mat_t res, const nmod_mat_t mat1, const nmod_mat_t mat2)
    void nmod_mat_concat_horizontal(nmod_mat_t res, const nmod_mat_t mat1, const nmod_mat_t mat2)
    void nmod_mat_print_pretty(nmod_mat_t mat)
    void nmod_mat_randtest(nmod_mat_t mat, flint_rand_t state)
    void nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state)
    int nmod_mat_randpermdiag(nmod_mat_t mat, mp_limb_t * diag, slong n, flint_rand_t state)
    void nmod_mat_randrank(nmod_mat_t mat, slong rank, flint_rand_t state)
    void nmod_mat_randops(nmod_mat_t mat, slong count, flint_rand_t state)
    void nmod_mat_randtril(nmod_mat_t mat, flint_rand_t state, int unit)
    void nmod_mat_randtriu(nmod_mat_t mat, flint_rand_t state, int unit)
    int nmod_mat_equal(nmod_mat_t mat1, nmod_mat_t mat2)
    int nmod_mat_is_zero_row(const nmod_mat_t mat, slong i)
    void nmod_mat_transpose(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_swap_rows(nmod_mat_t mat, slong * perm, slong r, slong s)
    void nmod_mat_swap_cols(nmod_mat_t mat, slong * perm, slong r, slong s)
    void nmod_mat_invert_rows(nmod_mat_t mat, slong * perm)
    void nmod_mat_invert_cols(nmod_mat_t mat, slong * perm)
    void nmod_mat_add(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_sub(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_neg(nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_scalar_mul(nmod_mat_t B, const nmod_mat_t A, mp_limb_t c)
    void nmod_mat_scalar_addmul_ui(nmod_mat_t dest, const nmod_mat_t X, const nmod_mat_t Y, const mp_limb_t b)
    void nmod_mat_scalar_mul_fmpz(nmod_mat_t res, const nmod_mat_t M, const fmpz_t c)
    void nmod_mat_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    void _nmod_mat_mul_classical_op(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, int op)
    void nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    # void _nmod_mat_mul_classical_threaded_pool_op(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, int op, thread_pool_handle * threads, slong num_threads)
    void _nmod_mat_mul_classical_threaded_op(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, int op)
    void nmod_mat_mul_classical_threaded(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    void nmod_mat_mul_strassen(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    void nmod_mat_addmul(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    void nmod_mat_submul(nmod_mat_t D, const nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
    void nmod_mat_mul_nmod_vec(mp_limb_t * c, const nmod_mat_t A, const mp_limb_t * b, slong blen)
    void nmod_mat_mul_nmod_vec_ptr(mp_limb_t * const * c, const nmod_mat_t A, const mp_limb_t * const * b, slong blen)
    void nmod_mat_nmod_vec_mul(mp_limb_t * c, const mp_limb_t * a, slong alen, const nmod_mat_t B)
    void nmod_mat_nmod_vec_mul_ptr(mp_limb_t * const * c, const mp_limb_t * const * a, slong alen, const nmod_mat_t B)
    void _nmod_mat_pow(nmod_mat_t dest, const nmod_mat_t mat, ulong pow)
    void nmod_mat_pow(nmod_mat_t dest, nmod_mat_t mat, ulong pow)
    mp_limb_t nmod_mat_trace(const nmod_mat_t mat)
    mp_limb_t nmod_mat_det_howell(const nmod_mat_t A)
    mp_limb_t nmod_mat_det(const nmod_mat_t A)
    slong nmod_mat_rank(const nmod_mat_t A)
    int nmod_mat_inv(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_solve_tril(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
    void nmod_mat_solve_tril_classical(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
    void nmod_mat_solve_tril_recursive(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit)
    void nmod_mat_solve_triu(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
    void nmod_mat_solve_triu_classical(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
    void nmod_mat_solve_triu_recursive(nmod_mat_t X, const nmod_mat_t U, const nmod_mat_t B, int unit)
    int nmod_mat_solve(nmod_mat_t X, nmod_mat_t A, nmod_mat_t B)
    int nmod_mat_can_solve_inner(slong * rank, slong * perm, slong * pivots, nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)
    int nmod_mat_can_solve(nmod_mat_t X, nmod_mat_t A, nmod_mat_t B)
    int nmod_mat_solve_vec(mp_limb_t * x, nmod_mat_t A, mp_limb_t * b)
    slong nmod_mat_lu(slong * P, nmod_mat_t A, int rank_check)
    slong nmod_mat_lu_classical(slong * P, nmod_mat_t A, int rank_check)
    slong nmod_mat_lu_classical_delayed(slong * P, nmod_mat_t A, int rank_check)
    slong nmod_mat_lu_recursive(slong * P, nmod_mat_t A, int rank_check)
    slong nmod_mat_rref(nmod_mat_t A)
    slong nmod_mat_reduce_row(nmod_mat_t A, slong * P, slong * L, slong n)
    slong nmod_mat_nullspace(nmod_mat_t X, const nmod_mat_t A)
    void nmod_mat_similarity(nmod_mat_t M, slong r, ulong d)
    void nmod_mat_charpoly_berkowitz(nmod_poly_t p, const nmod_mat_t M)
    void nmod_mat_charpoly_danilevsky(nmod_poly_t p, const nmod_mat_t M)
    void nmod_mat_charpoly(nmod_poly_t p, const nmod_mat_t M)
    void nmod_mat_minpoly(nmod_poly_t p, const nmod_mat_t M)
    void nmod_mat_strong_echelon_form(nmod_mat_t A)
    slong nmod_mat_howell_form(nmod_mat_t A)
    # not in doc
    int nmod_mat_is_square(nmod_mat_t mat)
