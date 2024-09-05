from flint.flintlib.types.flint cimport flint_bitcnt_t, flint_rand_t, fmpz_t, slong
from flint.flintlib.types.fmpz cimport fmpz_lll_t, fmpz_mat_t, gram_type, rep_type

# unknown type d_mat_t
# unknown type fmpz_gram_t
# unknown type mpf
# unknown type mpf_mat_t
# unknown type mpf_t


cdef extern from "flint/fmpz_lll.h":
    void fmpz_lll_context_init_default(fmpz_lll_t fl)
    void fmpz_lll_context_init(fmpz_lll_t fl, double delta, double eta, rep_type rt, gram_type gt)
    void fmpz_lll_randtest(fmpz_lll_t fl, flint_rand_t state)
    double fmpz_lll_heuristic_dot(const double * vec1, const double * vec2, slong len2, const fmpz_mat_t B, slong k, slong j, slong exp_adj)
    # int fmpz_lll_check_babai(int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
    # int fmpz_lll_check_babai_heuristic_d(int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
    # int fmpz_lll_check_babai_heuristic(int kappa, fmpz_mat_t B, fmpz_mat_t U, mpf_mat_t mu, mpf_mat_t r, mpf * s, mpf_mat_t appB, fmpz_gram_t A, int a, int zeros, int kappamax, int n, mpf_t tmp, mpf_t rtmp, flint_bitcnt_t prec, const fmpz_lll_t fl)
    # int fmpz_lll_advance_check_babai(int cur_kappa, int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
    # int fmpz_lll_advance_check_babai_heuristic_d(int cur_kappa, int kappa, fmpz_mat_t B, fmpz_mat_t U, d_mat_t mu, d_mat_t r, double * s, d_mat_t appB, int * expo, fmpz_gram_t A, int a, int zeros, int kappamax, int n, const fmpz_lll_t fl)
    int fmpz_lll_shift(const fmpz_mat_t B)
    int fmpz_lll_d(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    int fmpz_lll_d_heuristic(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    int fmpz_lll_mpf2(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_lll_t fl)
    int fmpz_lll_mpf(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    int fmpz_lll_wrapper(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    int fmpz_lll_d_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    int fmpz_lll_d_heuristic_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    int fmpz_lll_mpf2_with_removal(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_t gs_B, const fmpz_lll_t fl)
    int fmpz_lll_mpf_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    int fmpz_lll_wrapper_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    int fmpz_lll_d_with_removal_knapsack(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    int fmpz_lll_wrapper_with_removal_knapsack(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
    int fmpz_lll_with_removal_ulll(fmpz_mat_t FM, fmpz_mat_t UM, slong new_size, const fmpz_t gs_B, const fmpz_lll_t fl)
    int fmpz_lll_is_reduced_d(const fmpz_mat_t B, const fmpz_lll_t fl)
    int fmpz_lll_is_reduced_mpfr(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec)
    int fmpz_lll_is_reduced_d_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd)
    int fmpz_lll_is_reduced_mpfr_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd, flint_bitcnt_t prec)
    int fmpz_lll_is_reduced(const fmpz_mat_t B, const fmpz_lll_t fl, flint_bitcnt_t prec)
    int fmpz_lll_is_reduced_with_removal(const fmpz_mat_t B, const fmpz_lll_t fl, const fmpz_t gs_B, int newd, flint_bitcnt_t prec)
    void fmpz_lll_storjohann_ulll(fmpz_mat_t FM, slong new_size, const fmpz_lll_t fl)
    void fmpz_lll(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
    int fmpz_lll_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B, const fmpz_lll_t fl)
