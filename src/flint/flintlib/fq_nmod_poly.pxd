from flint.flintlib.fq_nmod cimport fq_nmod_struct, fq_nmod_ctx_t, fq_nmod_t, fq_nmod_poly_t, fq_nmod_mat_t
from flint.flintlib.flint cimport ulong, slong, flint_rand_t
from flint.flintlib.fmpz_mod_poly cimport fmpz_mod_poly_t
from flint.flintlib.nmod_poly cimport nmod_poly_t
from flint.flintlib.fmpz cimport fmpz_t

cdef extern from "flint/fq_nmod_poly.h":
    # Parsed from here **********************************************/
    void fq_nmod_poly_init(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_init2(fq_nmod_poly_t poly, slong alloc, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_realloc(fq_nmod_poly_t poly, slong alloc, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_fit_length(fq_nmod_poly_t poly, slong len, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_set_length(fq_nmod_poly_t poly, slong newlen, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_clear(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_normalise(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_normalise2(const fq_nmod_struct * poly, slong * length, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_truncate(fq_nmod_poly_t poly, slong newlen, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_set_trunc(fq_nmod_poly_t poly1, fq_nmod_poly_t poly2, slong newlen, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_reverse(fq_nmod_struct * output, const fq_nmod_struct * input, slong len, slong m, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_reverse(fq_nmod_poly_t output, const fq_nmod_poly_t input, slong m, const fq_nmod_ctx_t ctx)
    slong fq_nmod_poly_degree(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    slong fq_nmod_poly_length(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    fq_nmod_struct * fq_nmod_poly_lead(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_randtest(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_randtest_not_zero(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_randtest_monic(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_randtest_irreducible(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_set(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_set(fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_set_fq_nmod(fq_nmod_poly_t poly, const fq_nmod_t c, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_set_fmpz_mod_poly(fq_nmod_poly_t rop, const fmpz_mod_poly_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_set_nmod_poly(fq_nmod_poly_t rop, const nmod_poly_t op, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_swap(fq_nmod_poly_t op1, fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_zero(fq_nmod_struct * rop, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_zero(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_one(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_gen(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_make_monic(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_make_monic(fq_nmod_struct * rop, const fq_nmod_struct * op, slong length, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_get_coeff(fq_nmod_t x, const fq_nmod_poly_t poly, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_set_coeff(fq_nmod_poly_t poly, slong n, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_set_coeff_fmpz(fq_nmod_poly_t poly, slong n, const fmpz_t x, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_equal(const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_equal_trunc(const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong n, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_is_zero(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_is_one(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_is_gen(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_is_unit(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_equal_fq_nmod(const fq_nmod_poly_t poly, const fq_nmod_t c, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_add(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_add(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_add_si(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, slong c, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_add_series(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong n, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_sub(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_sub(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_sub_series(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong n, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_neg(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_neg(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_scalar_mul_fq_nmod(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_scalar_mul_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_scalar_addmul_fq_nmod(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_scalar_addmul_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_scalar_submul_fq_nmod(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_scalar_submul_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_scalar_div_fq(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_scalar_div_fq(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mul_classical(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mul_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mul_reorder(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mul_reorder(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mul_univariate(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mul_univariate(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mul_KS(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mul_KS(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mul(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mul(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mullow_classical(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mullow_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mullow_univariate(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mullow_univariate(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mullow_KS(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mullow_KS(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mullow(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mullow(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mulhigh_classical(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, slong start, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mulhigh_classical(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong start, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mulhigh(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, slong start, fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mulhigh(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong start, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mulmod(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, const fq_nmod_struct * f, slong lenf, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mulmod(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_mulmod_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_mulmod_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_sqr_classical(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_sqr_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_sqr_KS(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_sqr_KS(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_sqr(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_sqr(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_pow(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, ulong e, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_pow(fq_nmod_poly_t rop, const fq_nmod_poly_t op, ulong e, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_powmod_ui_binexp(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, const fq_nmod_struct * f, slong lenf, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_powmod_ui_binexp(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_powmod_ui_binexp_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_powmod_ui_binexp_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_powmod_fmpz_binexp(fq_nmod_struct * res, const fq_nmod_struct * poly, const fmpz_t e, const fq_nmod_struct * f, slong lenf, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_powmod_fmpz_binexp(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fmpz_t e, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_powmod_fmpz_binexp_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly, const fmpz_t e, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_powmod_fmpz_binexp_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fmpz_t e, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_powmod_fmpz_sliding_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly, const fmpz_t e, ulong k, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_powmod_fmpz_sliding_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fmpz_t e, ulong k, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_powmod_x_fmpz_preinv(fq_nmod_struct * res, const fmpz_t e, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_powmod_x_fmpz_preinv(fq_nmod_poly_t res, const fmpz_t e, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_pow_trunc_binexp(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_pow_trunc_binexp(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_pow_trunc(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, slong trunc, const fq_nmod_ctx_t mod)
    void fq_nmod_poly_pow_trunc(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_shift_left(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_shift_left(fq_nmod_poly_t rop, const fq_nmod_poly_t op, slong n, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_shift_right(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_shift_right(fq_nmod_poly_t rop, const fq_nmod_poly_t op, slong n, const fq_nmod_ctx_t ctx)
    slong _fq_nmod_poly_hamming_weight(const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
    slong fq_nmod_poly_hamming_weight(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_divrem(fq_nmod_struct * Q, fq_nmod_struct * R, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_divrem(fq_nmod_poly_t Q, fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_divrem_f(fq_nmod_t f, fq_nmod_poly_t Q, fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_rem(fq_nmod_struct * R, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_rem(fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_div(fq_nmod_struct * Q, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_div(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_div_newton_n_preinv(fq_nmod_struct * Q, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_struct * Binv, slong lenBinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_div_newton_n_preinv(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_poly_t Binv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_divrem_newton_n_preinv(fq_nmod_struct * Q, fq_nmod_struct * R, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_struct * Binv, slong lenBinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_divrem_newton_n_preinv(fq_nmod_poly_t Q, fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_poly_t Binv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_inv_series_newton(fq_nmod_struct * Qinv, const fq_nmod_struct * Q, slong n, const fq_nmod_t cinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_inv_series_newton(fq_nmod_poly_t Qinv, const fq_nmod_poly_t Q, slong n, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_inv_series(fq_nmod_struct * Qinv, const fq_nmod_struct * Q, slong n, const fq_nmod_t cinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_inv_series(fq_nmod_poly_t Qinv, const fq_nmod_poly_t Q, slong n, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_div_series(fq_nmod_struct * Q, const fq_nmod_struct * A, slong Alen, const fq_nmod_struct * B, slong Blen, slong n, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_div_series(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, slong n, fq_nmod_ctx_t ctx)
    void fq_nmod_poly_gcd(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    slong _fq_nmod_poly_gcd(fq_nmod_struct * G, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_ctx_t ctx)
    slong _fq_nmod_poly_gcd_euclidean_f(fq_nmod_t f, fq_nmod_struct * G, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_gcd_euclidean_f(fq_nmod_t f, fq_nmod_poly_t G, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    slong _fq_nmod_poly_xgcd(fq_nmod_struct * G, fq_nmod_struct * S, fq_nmod_struct * T, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_xgcd(fq_nmod_poly_t G, fq_nmod_poly_t S, fq_nmod_poly_t T, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    slong _fq_nmod_poly_xgcd_euclidean_f(fq_nmod_t f, fq_nmod_struct * G, fq_nmod_struct * S, fq_nmod_struct * T, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_xgcd_euclidean_f(fq_nmod_t f, fq_nmod_poly_t G, fq_nmod_poly_t S, fq_nmod_poly_t T, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    int _fq_nmod_poly_divides(fq_nmod_struct * Q, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_divides(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_derivative(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_derivative(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_invsqrt_series(fq_nmod_struct * g, const fq_nmod_struct * h, slong n, fq_nmod_ctx_t mod)
    void fq_nmod_poly_invsqrt_series(fq_nmod_poly_t g, const fq_nmod_poly_t h, slong n, fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_sqrt_series(fq_nmod_struct * g, const fq_nmod_struct * h, slong n, fq_nmod_ctx_t ctx)
    void fq_nmod_poly_sqrt_series(fq_nmod_poly_t g, const fq_nmod_poly_t h, slong n, fq_nmod_ctx_t ctx)
    int _fq_nmod_poly_sqrt(fq_nmod_struct * s, const fq_nmod_struct * p, slong n, fq_nmod_ctx_t mod)
    int fq_nmod_poly_sqrt(fq_nmod_poly_t s, const fq_nmod_poly_t p, fq_nmod_ctx_t mod)
    void _fq_nmod_poly_evaluate_fq_nmod(fq_nmod_t rop, const fq_nmod_struct * op, slong len, const fq_nmod_t a, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_evaluate_fq_nmod(fq_nmod_t rop, const fq_nmod_poly_t f, const fq_nmod_t a, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_compose(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_compose(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_compose_mod_horner(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_compose_mod_horner(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_compose_mod_horner_preinv(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_struct * hinv, slong lenhiv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_compose_mod_horner_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_poly_t hinv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_compose_mod_brent_kung(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_compose_mod_brent_kung(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_compose_mod_brent_kung_preinv(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_struct * hinv, slong lenhiv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_compose_mod_brent_kung_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_poly_t hinv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_compose_mod(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_compose_mod(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_compose_mod_preinv(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_struct * hinv, slong lenhiv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_compose_mod_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t h, const fq_nmod_poly_t hinv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_reduce_matrix_mod_poly (fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_precompute_matrix (fq_nmod_mat_t A, const fq_nmod_struct * f, const fq_nmod_struct * g, slong leng, const fq_nmod_struct * ginv, slong lenginv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_precompute_matrix (fq_nmod_mat_t A, const fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_poly_t ginv, const fq_nmod_ctx_t ctx)
    void _fq_nmod_poly_compose_mod_brent_kung_precomp_preinv(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_mat_t A, const fq_nmod_struct * h, slong lenh, const fq_nmod_struct * hinv, slong lenhinv, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_compose_mod_brent_kung_precomp_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t f, const fq_nmod_mat_t A, const fq_nmod_poly_t h, const fq_nmod_poly_t hinv, const fq_nmod_ctx_t ctx)
    # int _fq_nmod_poly_fprint_pretty(FILE * file, const fq_nmod_struct * poly, slong len, const char * x, const fq_nmod_ctx_t ctx)
    # int fq_nmod_poly_fprint_pretty(FILE * file, const fq_nmod_poly_t poly, const char * x, const fq_nmod_ctx_t ctx)
    int _fq_nmod_poly_print_pretty(const fq_nmod_struct * poly, slong len, const char * x, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_print_pretty(const fq_nmod_poly_t poly, const char * x, const fq_nmod_ctx_t ctx)
    # int _fq_nmod_poly_fprint(FILE * file, const fq_nmod_struct * poly, slong len, const fq_nmod_ctx_t ctx)
    # int fq_nmod_poly_fprint(FILE * file, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    int _fq_nmod_poly_print(const fq_nmod_struct * poly, slong len, const fq_nmod_ctx_t ctx)
    int fq_nmod_poly_print(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    char * _fq_nmod_poly_get_str(const fq_nmod_struct * poly, slong len, const fq_nmod_ctx_t ctx)
    char * fq_nmod_poly_get_str(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    char * _fq_nmod_poly_get_str_pretty(const fq_nmod_struct * poly, slong len, const char * x, const fq_nmod_ctx_t ctx)
    char * fq_nmod_poly_get_str_pretty(const fq_nmod_poly_t poly, const char * x, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_inflate(fq_nmod_poly_t result, const fq_nmod_poly_t input, ulong inflation, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_deflate(fq_nmod_poly_t result, const fq_nmod_poly_t input, ulong deflation, const fq_nmod_ctx_t ctx)
    ulong fq_nmod_poly_deflation(const fq_nmod_poly_t input, const fq_nmod_ctx_t ctx)