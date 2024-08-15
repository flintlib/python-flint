from flint.flintlib.flint cimport flint_rand_t, slong
from flint.flintlib.fq cimport fq_ctx_t, fq_struct
from flint.flintlib.nmod_mat cimport nmod_mat_t
from flint.flintlib.fmpz_mat cimport fmpz_mat_t
from flint.flintlib.fmpz_mod_mat cimport fmpz_mod_mat_t
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.fq_default cimport fq_default_t, fq_default_ctx_t
from flint.flintlib.fq_default_poly cimport fq_default_poly_t
from flint.flintlib.fq_mat cimport fq_mat_t
from flint.flintlib.fq_nmod_mat cimport fq_nmod_mat_t
from flint.flintlib.fq_zech_mat cimport fq_zech_mat_t
from flintlib.nmod_mat cimport nmod_mat_t
from flintlib.fmpz_mod_mat cimport fmpz_mod_mat_t

cdef extern from "flint/fq_default_mat.h":
    # Type definitions **********************************************/
    ctypedef union fq_default_mat_struct:
        fq_mat_t fq
        fq_nmod_mat_t fq_nmod
        fq_zech_mat_t fq_zech
        nmod_mat_t nmod
        fmpz_mod_mat_t fmpz_mod
    ctypedef fq_default_mat_struct fq_default_mat_t[1]

    # Parsed from here **********************************************/
    void fq_default_mat_init(fq_default_mat_t mat, slong rows, slong cols, const fq_default_ctx_t ctx)
    void fq_default_mat_init_set(fq_default_mat_t mat, const fq_default_mat_t src, const fq_default_ctx_t ctx)
    void fq_default_mat_clear(fq_default_mat_t mat, const fq_default_ctx_t ctx)
    void fq_default_mat_set(fq_default_mat_t mat, const fq_default_mat_t src, const fq_default_ctx_t ctx)
    void fq_default_mat_entry(fq_default_t val, const fq_default_mat_t mat, slong i, slong j, const fq_default_ctx_t ctx)
    void fq_default_mat_entry_set(fq_default_mat_t mat, slong i, slong j, const fq_default_t x, const fq_default_ctx_t ctx)
    void fq_default_mat_entry_set_fmpz(fq_default_mat_t mat, slong i, slong j, const fmpz_t x, const fq_default_ctx_t ctx)
    slong fq_default_mat_nrows(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    slong fq_default_mat_ncols(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    void fq_default_mat_swap(fq_default_mat_t mat1, fq_default_mat_t mat2, const fq_default_ctx_t ctx)
    void fq_default_mat_zero(fq_default_mat_t mat, const fq_default_ctx_t ctx)
    void fq_default_mat_one(fq_default_mat_t mat, const fq_default_ctx_t ctx)
    void fq_default_mat_swap_rows(fq_default_mat_t mat, slong * perm, slong r, slong s, const fq_default_ctx_t ctx)
    void fq_default_mat_swap_cols(fq_default_mat_t mat, slong * perm, slong r, slong s, const fq_default_ctx_t ctx)
    void fq_default_mat_invert_rows(fq_default_mat_t mat, slong * perm, const fq_default_ctx_t ctx)
    void fq_default_mat_invert_cols(fq_default_mat_t mat, slong * perm, const fq_default_ctx_t ctx)
    void fq_default_mat_set_nmod_mat(fq_default_mat_t mat1, const nmod_mat_t mat2, const fq_default_ctx_t ctx)
    void fq_default_mat_set_fmpz_mod_mat(fq_default_mat_t mat1, const fmpz_mod_mat_t mat2, const fq_default_ctx_t ctx)
    void fq_default_mat_set_fmpz_mat(fq_default_mat_t mat1, const fmpz_mat_t mat2, const fq_default_ctx_t ctx)
    void fq_default_mat_concat_vertical(fq_default_mat_t res, const fq_default_mat_t mat1, const fq_default_mat_t mat2, const fq_default_ctx_t ctx)
    void fq_default_mat_concat_horizontal(fq_default_mat_t res, const fq_default_mat_t mat1, const fq_default_mat_t mat2, const fq_default_ctx_t ctx)
    int fq_default_mat_print_pretty(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    # int fq_default_mat_fprint_pretty(FILE * file, const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    int fq_default_mat_print(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    # int fq_default_mat_fprint(FILE * file, const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    void fq_default_mat_window_init(fq_default_mat_t window, const fq_default_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fq_default_ctx_t ctx)
    void fq_default_mat_window_clear(fq_default_mat_t window, const fq_default_ctx_t ctx)
    void fq_default_mat_randtest(fq_default_mat_t mat, flint_rand_t state, const fq_default_ctx_t ctx)
    int fq_default_mat_randpermdiag(fq_mat_t mat, flint_rand_t state, fq_struct * diag, slong n, const fq_ctx_t ctx)
    void fq_default_mat_randrank(fq_default_mat_t mat, flint_rand_t state, slong rank, const fq_default_ctx_t ctx)
    void fq_default_mat_randops(fq_default_mat_t mat, flint_rand_t state, slong count, const fq_default_ctx_t ctx)
    void fq_default_mat_randtril(fq_default_mat_t mat, flint_rand_t state, int unit, const fq_default_ctx_t ctx)
    void fq_default_mat_randtriu(fq_default_mat_t mat, flint_rand_t state, int unit, const fq_default_ctx_t ctx)
    int fq_default_mat_equal(const fq_default_mat_t mat1, const fq_default_mat_t mat2, const fq_default_ctx_t ctx)
    int fq_default_mat_is_zero(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    int fq_default_mat_is_one(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    int fq_default_mat_is_empty(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    int fq_default_mat_is_square(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
    void fq_default_mat_add(fq_default_mat_t C, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)
    void fq_default_mat_sub(fq_default_mat_t C, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)
    void fq_default_mat_neg(fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)
    void fq_default_mat_mul(fq_default_mat_t C, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)
    void fq_default_mat_submul(fq_default_mat_t D, const fq_default_mat_t C, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)
    int fq_default_mat_inv(fq_default_mat_t B, fq_default_mat_t A, const fq_default_ctx_t ctx)
    slong fq_default_mat_lu(slong * P, fq_default_mat_t A, int rank_check, const fq_default_ctx_t ctx)
    slong fq_default_mat_rref(fq_default_mat_t B, const fq_default_mat_t A, const fq_default_ctx_t ctx)
    void fq_default_mat_solve_tril(fq_default_mat_t X, const fq_default_mat_t L, const fq_default_mat_t B, int unit, const fq_default_ctx_t ctx)
    void fq_default_mat_solve_triu(fq_default_mat_t X, const fq_default_mat_t U, const fq_default_mat_t B, int unit, const fq_default_ctx_t ctx)
    int fq_default_mat_solve(fq_default_mat_t X, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)
    int fq_default_mat_can_solve(fq_default_mat_t X, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)
    void fq_default_mat_similarity(fq_default_mat_t M, slong r, fq_default_t d, const fq_default_ctx_t ctx)
    void fq_default_mat_charpoly(fq_default_poly_t p, const fq_default_mat_t M, const fq_default_ctx_t ctx)
    void fq_default_mat_minpoly(fq_default_poly_t p, const fq_default_mat_t M, const fq_default_ctx_t ctx)