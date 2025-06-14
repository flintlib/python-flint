from flint.flintlib.types.flint cimport flint_rand_t, slong
from flint.flintlib.types.fmpz_mod cimport fmpz_mod_mat_t
from flint.flintlib.types.fq_nmod cimport fq_nmod_ctx_t, fq_nmod_mat_t, fq_nmod_poly_t, fq_nmod_struct, fq_nmod_t
from flint.flintlib.types.nmod cimport nmod_mat_t

# unknown type FILE


cdef extern from "flint/fq_nmod_mat.h":
    void fq_nmod_mat_init(fq_nmod_mat_t mat, slong rows, slong cols, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_init_set(fq_nmod_mat_t mat, const fq_nmod_mat_t src, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_clear(fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_set(fq_nmod_mat_t mat, const fq_nmod_mat_t src, const fq_nmod_ctx_t ctx)
    fq_nmod_struct * fq_nmod_mat_entry(const fq_nmod_mat_t mat, slong i, slong j)
    void fq_nmod_mat_entry_set(fq_nmod_mat_t mat, slong i, slong j, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    slong fq_nmod_mat_nrows(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    slong fq_nmod_mat_ncols(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_swap(fq_nmod_mat_t mat1, fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_swap_entrywise(fq_nmod_mat_t mat1, fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_zero(fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_one(fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_swap_rows(fq_nmod_mat_t mat, slong * perm, slong r, slong s, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_swap_cols(fq_nmod_mat_t mat, slong * perm, slong r, slong s, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_invert_rows(fq_nmod_mat_t mat, slong * perm, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_invert_cols(fq_nmod_mat_t mat, slong * perm, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_set_nmod_mat(fq_nmod_mat_t mat1, const nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_set_fmpz_mod_mat(fq_nmod_mat_t mat1, const fmpz_mod_mat_t mat2, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_concat_vertical(fq_nmod_mat_t res, const fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_concat_horizontal(fq_nmod_mat_t res, const fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_print_pretty(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # int fq_nmod_mat_fprint_pretty(FILE * file, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_print(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    # int fq_nmod_mat_fprint(FILE * file, const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_window_init(fq_nmod_mat_t window, const fq_nmod_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_window_clear(fq_nmod_mat_t window, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_randtest(fq_nmod_mat_t mat, flint_rand_t state, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_randpermdiag(fq_nmod_mat_t mat, flint_rand_t state, fq_nmod_struct * diag, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_randrank(fq_nmod_mat_t mat, flint_rand_t state, slong rank, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_randops(fq_nmod_mat_t mat, flint_rand_t state, slong count, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_randtril(fq_nmod_mat_t mat, flint_rand_t state, int unit, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_randtriu(fq_nmod_mat_t mat, flint_rand_t state, int unit, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_equal(const fq_nmod_mat_t mat1, const fq_nmod_mat_t mat2, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_is_zero(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_is_one(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_is_empty(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_is_square(const fq_nmod_mat_t mat, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_transpose(fq_nmod_mat_t B, const fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_add(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B,  const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_sub(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_neg(fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_mul(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B,  const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_mul_classical(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_mul_KS(fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_submul(fq_nmod_mat_t D, const fq_nmod_mat_t C, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_mul_vec(fq_nmod_struct * c, const fq_nmod_mat_t A, const fq_nmod_struct * b, slong blen, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_mul_vec_ptr(fq_nmod_struct * const * c, const fq_nmod_mat_t A, const fq_nmod_struct * const * b, slong blen, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_vec_mul(fq_nmod_struct * c, const fq_nmod_struct * a, slong alen, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_vec_mul_ptr(fq_nmod_struct * const * c, const fq_nmod_struct * const * a, slong alen, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_inv(fq_nmod_mat_t B, fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)
    slong fq_nmod_mat_lu(slong * P, fq_nmod_mat_t A, int rank_check, const fq_nmod_ctx_t ctx)
    slong fq_nmod_mat_rref(fq_nmod_mat_t B, const fq_nmod_mat_t A, const fq_nmod_ctx_t ctx)
    slong fq_nmod_mat_reduce_row(fq_nmod_mat_t A, slong * P, slong * L, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_solve_tril(fq_nmod_mat_t X, const fq_nmod_mat_t L, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_solve_triu(fq_nmod_mat_t X, const fq_nmod_mat_t U, const fq_nmod_mat_t B, int unit, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_solve(fq_nmod_mat_t X, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    int fq_nmod_mat_can_solve(fq_nmod_mat_t X, const fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_similarity(fq_nmod_mat_t M, slong r, fq_nmod_t d, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_charpoly(fq_nmod_poly_t p, const fq_nmod_mat_t M, const fq_nmod_ctx_t ctx)
    void fq_nmod_mat_minpoly(fq_nmod_poly_t p, const fq_nmod_mat_t M, const fq_nmod_ctx_t ctx)
