from flint.flintlib.types.flint cimport flint_bitcnt_t, flint_rand_t, fmpz_struct, fmpz_t, slong, ulong
from flint.flintlib.types.fmpq cimport fmpq_t
from flint.flintlib.types.gr cimport gr_ctx_t, gr_mpoly_t, gr_ptr, gr_srcptr, gr_stream_t, truth_t
from flint.flintlib.types.mpoly cimport mpoly_ctx_t



cdef extern from "flint/gr_mpoly.h":
    void gr_mpoly_init(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_init3(gr_mpoly_t A, slong alloc, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_init2(gr_mpoly_t A, slong alloc, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_clear(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_swap(gr_mpoly_t A, gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set(gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_zero(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    truth_t gr_mpoly_is_zero(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_gen(gr_mpoly_t A, slong var, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    truth_t gr_mpoly_is_gen(const gr_mpoly_t A, slong var, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    truth_t gr_mpoly_equal(const gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_randtest_bits(gr_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_write_pretty(gr_stream_t out, const gr_mpoly_t A, const char ** x, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_print_pretty(const gr_mpoly_t A, const char ** x, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_get_coeff_scalar_fmpz(gr_ptr c, const gr_mpoly_t A, const fmpz_struct * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_get_coeff_scalar_ui(gr_ptr c, const gr_mpoly_t A, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz_struct * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_ui_fmpz(gr_mpoly_t A, ulong c, const fmpz_struct * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_si_fmpz(gr_mpoly_t A, slong c, const fmpz_struct * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_fmpz_fmpz(gr_mpoly_t A, const fmpz_t c, const fmpz_struct * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_fmpq_fmpz(gr_mpoly_t A, const fmpq_t c, const fmpz_struct * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_scalar_ui(gr_mpoly_t poly, gr_srcptr c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_ui_ui(gr_mpoly_t A, ulong c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_si_ui(gr_mpoly_t A, slong c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_fmpz_ui(gr_mpoly_t A, const fmpz_t c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_set_coeff_fmpq_ui(gr_mpoly_t A, const fmpq_t c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_neg(gr_mpoly_t A, const gr_mpoly_t B, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_add(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_sub(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_johnson(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_monomial(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_scalar(gr_mpoly_t A, const gr_mpoly_t B, gr_srcptr c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_si(gr_mpoly_t A, const gr_mpoly_t B, slong c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_ui(gr_mpoly_t A, const gr_mpoly_t B, ulong c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_fmpz(gr_mpoly_t A, const gr_mpoly_t B, const fmpz_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_mul_fmpq(gr_mpoly_t A, const gr_mpoly_t B, const fmpq_t c, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void _gr_mpoly_fit_length(gr_ptr * coeffs, slong * coeffs_alloc, ulong ** exps, slong * exps_alloc, slong N, slong length, gr_ctx_t cctx)
    void gr_mpoly_fit_length(gr_mpoly_t A, slong len, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_fit_bits(gr_mpoly_t A, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_fit_length_fit_bits(gr_mpoly_t A, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_fit_length_reset_bits(gr_mpoly_t A, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void _gr_mpoly_set_length(gr_mpoly_t A, slong newlen, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void _gr_mpoly_push_exp_ui(gr_mpoly_t A, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_push_term_scalar_ui(gr_mpoly_t A, gr_srcptr c, const ulong * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void _gr_mpoly_push_exp_fmpz(gr_mpoly_t A, const fmpz_struct * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_push_term_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz_struct * exp, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_sort_terms(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    int gr_mpoly_combine_like_terms(gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    truth_t gr_mpoly_is_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)
    void gr_mpoly_assert_canonical(const gr_mpoly_t A, const mpoly_ctx_t mctx, gr_ctx_t cctx)