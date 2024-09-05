from flint.flintlib.types.flint cimport flint_bitcnt_t, flint_rand_t, fmpz_struct, fmpz_t, slong, ulong
from flint.flintlib.types.fmpq cimport fmpq_poly_struct, fmpq_poly_t, fmpq_t
from flint.flintlib.types.fmpz cimport fmpz_poly_t, fmpz_preinvn_t
from flint.flintlib.types.nmod cimport nmod_poly_t

# unknown type FILE
# unknown type fmpq_poly_powers_precomp_t


cdef extern from "flint/fmpq_poly.h":
    void fmpq_poly_init(fmpq_poly_t poly)
    void fmpq_poly_init2(fmpq_poly_t poly, slong alloc)
    void fmpq_poly_realloc(fmpq_poly_t poly, slong alloc)
    void fmpq_poly_fit_length(fmpq_poly_t poly, slong len)
    void _fmpq_poly_set_length(fmpq_poly_t poly, slong len)
    void fmpq_poly_clear(fmpq_poly_t poly)
    void _fmpq_poly_normalise(fmpq_poly_t poly)
    void _fmpq_poly_canonicalise(fmpz_struct * poly, fmpz_t den, slong len)
    void fmpq_poly_canonicalise(fmpq_poly_t poly)
    int _fmpq_poly_is_canonical(const fmpz_struct * poly, const fmpz_t den, slong len)
    int fmpq_poly_is_canonical(const fmpq_poly_t poly)
    slong fmpq_poly_degree(const fmpq_poly_t poly)
    slong fmpq_poly_length(const fmpq_poly_t poly)
    fmpz_struct * fmpq_poly_numref(fmpq_poly_t poly)
    fmpz_t fmpq_poly_denref(fmpq_poly_t poly)
    void fmpq_poly_get_numerator(fmpz_poly_t res, const fmpq_poly_t poly)
    void fmpq_poly_get_denominator(fmpz_t den, const fmpq_poly_t poly)
    void fmpq_poly_randtest(fmpq_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)
    void fmpq_poly_randtest_unsigned(fmpq_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)
    void fmpq_poly_randtest_not_zero(fmpq_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)
    void fmpq_poly_set(fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void fmpq_poly_set_si(fmpq_poly_t poly, slong x)
    void fmpq_poly_set_ui(fmpq_poly_t poly, ulong x)
    void fmpq_poly_set_fmpz(fmpq_poly_t poly, const fmpz_t x)
    void fmpq_poly_set_fmpq(fmpq_poly_t poly, const fmpq_t x)
    void fmpq_poly_set_fmpz_poly(fmpq_poly_t rop, const fmpz_poly_t op)
    void fmpq_poly_set_nmod_poly(fmpq_poly_t rop, const nmod_poly_t op)
    void fmpq_poly_get_nmod_poly(nmod_poly_t rop, const fmpq_poly_t op)
    void fmpq_poly_get_nmod_poly_den(nmod_poly_t rop, const fmpq_poly_t op, int den)
    int _fmpq_poly_set_str(fmpz_struct * poly, fmpz_t den, const char * str, slong len)
    int fmpq_poly_set_str(fmpq_poly_t poly, const char * str)
    char * fmpq_poly_get_str(const fmpq_poly_t poly)
    char * fmpq_poly_get_str_pretty(const fmpq_poly_t poly, const char * var)
    void fmpq_poly_zero(fmpq_poly_t poly)
    void fmpq_poly_one(fmpq_poly_t poly)
    void fmpq_poly_neg(fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void fmpq_poly_inv(fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void fmpq_poly_swap(fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_truncate(fmpq_poly_t poly, slong n)
    void fmpq_poly_set_trunc(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void fmpq_poly_get_slice(fmpq_poly_t rop, const fmpq_poly_t op, slong i, slong j)
    void fmpq_poly_reverse(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void fmpq_poly_get_coeff_fmpz(fmpz_t x, const fmpq_poly_t poly, slong n)
    void fmpq_poly_get_coeff_fmpq(fmpq_t x, const fmpq_poly_t poly, slong n)
    void fmpq_poly_set_coeff_si(fmpq_poly_t poly, slong n, slong x)
    void fmpq_poly_set_coeff_ui(fmpq_poly_t poly, slong n, ulong x)
    void fmpq_poly_set_coeff_fmpz(fmpq_poly_t poly, slong n, const fmpz_t x)
    void fmpq_poly_set_coeff_fmpq(fmpq_poly_t poly, slong n, const fmpq_t x)
    int fmpq_poly_equal(const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    int _fmpq_poly_equal_trunc(const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, slong n)
    int fmpq_poly_equal_trunc(const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    int _fmpq_poly_cmp(const fmpz_struct * lpoly, const fmpz_t lden, const fmpz_struct * rpoly, const fmpz_t rden, slong len)
    int fmpq_poly_cmp(const fmpq_poly_t left, const fmpq_poly_t right)
    int fmpq_poly_is_one(const fmpq_poly_t poly)
    int fmpq_poly_is_zero(const fmpq_poly_t poly)
    int fmpq_poly_is_gen(const fmpq_poly_t poly)
    void _fmpq_poly_add(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2)
    void _fmpq_poly_add_can(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, int can)
    void fmpq_poly_add(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void fmpq_poly_add_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, int can)
    void _fmpq_poly_add_series(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, slong n)
    void _fmpq_poly_add_series_can(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, slong n, int can)
    void fmpq_poly_add_series(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    void fmpq_poly_add_series_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n, int can)
    void _fmpq_poly_sub(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2)
    void _fmpq_poly_sub_can(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, int can)
    void fmpq_poly_sub(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void fmpq_poly_sub_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, int can)
    void _fmpq_poly_sub_series(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, slong n)
    void _fmpq_poly_sub_series_can(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, slong n, int can)
    void fmpq_poly_sub_series(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    void fmpq_poly_sub_series_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n, int can)
    void _fmpq_poly_scalar_mul_si(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, slong c)
    void _fmpq_poly_scalar_mul_ui(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, ulong c)
    void _fmpq_poly_scalar_mul_fmpz(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, const fmpz_t c)
    void _fmpq_poly_scalar_mul_fmpq(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, const fmpz_t r, const fmpz_t s)
    void fmpq_poly_scalar_mul_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
    void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c)
    void fmpq_poly_scalar_mul_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
    void fmpq_poly_scalar_mul_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
    void _fmpq_poly_scalar_div_fmpz(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, const fmpz_t c)
    void _fmpq_poly_scalar_div_si(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, slong c)
    void _fmpq_poly_scalar_div_ui(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, ulong c)
    void _fmpq_poly_scalar_div_fmpq(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, const fmpz_t r, const fmpz_t s)
    void fmpq_poly_scalar_div_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c)
    void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
    void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
    void fmpq_poly_scalar_div_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
    void _fmpq_poly_mul(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2)
    void fmpq_poly_mul(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void _fmpq_poly_mullow(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, slong n)
    void fmpq_poly_mullow(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    void fmpq_poly_addmul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2)
    void fmpq_poly_submul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2)
    void _fmpq_poly_pow(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, ulong e)
    void fmpq_poly_pow(fmpq_poly_t res, const fmpq_poly_t poly, ulong e)
    void _fmpq_poly_pow_trunc(fmpz_struct * res, fmpz_t rden, const fmpz_struct * f, const fmpz_t fden, slong flen, ulong exp, slong len)
    void fmpq_poly_pow_trunc(fmpq_poly_t res, const fmpq_poly_t poly, ulong e, slong n)
    void fmpq_poly_shift_left(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void fmpq_poly_shift_right(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void _fmpq_poly_divrem(fmpz_struct * Q, fmpz_t q, fmpz_struct * R, fmpz_t r, const fmpz_struct * A, const fmpz_t a, slong lenA, const fmpz_struct * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)
    void fmpq_poly_divrem(fmpq_poly_t Q, fmpq_poly_t R, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void _fmpq_poly_div(fmpz_struct * Q, fmpz_t q, const fmpz_struct * A, const fmpz_t a, slong lenA, const fmpz_struct * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)
    void fmpq_poly_div(fmpq_poly_t Q, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void _fmpq_poly_rem(fmpz_struct * R, fmpz_t r, const fmpz_struct * A, const fmpz_t a, slong lenA, const fmpz_struct * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)
    void fmpq_poly_rem(fmpq_poly_t R, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    fmpq_poly_struct * _fmpq_poly_powers_precompute(const fmpz_struct * B, const fmpz_t denB, slong len)
    # void fmpq_poly_powers_precompute(fmpq_poly_powers_precomp_t pinv, fmpq_poly_t poly)
    void _fmpq_poly_powers_clear(fmpq_poly_struct * powers, slong len)
    # void fmpq_poly_powers_clear(fmpq_poly_powers_precomp_t pinv)
    void _fmpq_poly_rem_powers_precomp(fmpz_struct * A, fmpz_t denA, slong m, const fmpz_struct * B, const fmpz_t denB, slong n, fmpq_poly_struct * const powers)
    # void fmpq_poly_rem_powers_precomp(fmpq_poly_t R, const fmpq_poly_t A, const fmpq_poly_t B, const fmpq_poly_powers_precomp_t B_inv)
    int _fmpq_poly_divides(fmpz_struct * qpoly, fmpz_t qden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2)
    int fmpq_poly_divides(fmpq_poly_t q, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    slong fmpq_poly_remove(fmpq_poly_t q, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void _fmpq_poly_inv_series_newton(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, slong n)
    void fmpq_poly_inv_series_newton(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void _fmpq_poly_inv_series(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong den_len, slong n)
    void fmpq_poly_inv_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void _fmpq_poly_div_series(fmpz_struct * Q, fmpz_t denQ, const fmpz_struct * A, const fmpz_t denA, slong lenA, const fmpz_struct * B, const fmpz_t denB, slong lenB, slong n)
    void fmpq_poly_div_series(fmpq_poly_t Q, const fmpq_poly_t A, const fmpq_poly_t B, slong n)
    void _fmpq_poly_gcd(fmpz_struct * G, fmpz_t denG, const fmpz_struct * A, slong lenA, const fmpz_struct * B, slong lenB)
    void fmpq_poly_gcd(fmpq_poly_t G, const fmpq_poly_t A, const fmpq_poly_t B)
    void _fmpq_poly_xgcd(fmpz_struct * G, fmpz_t denG, fmpz_struct * S, fmpz_t denS, fmpz_struct * T, fmpz_t denT, const fmpz_struct * A, const fmpz_t denA, slong lenA, const fmpz_struct * B, const fmpz_t denB, slong lenB)
    void fmpq_poly_xgcd(fmpq_poly_t G, fmpq_poly_t S, fmpq_poly_t T, const fmpq_poly_t A, const fmpq_poly_t B)
    void _fmpq_poly_lcm(fmpz_struct * L, fmpz_t denL, const fmpz_struct * A, slong lenA, const fmpz_struct * B, slong lenB)
    void fmpq_poly_lcm(fmpq_poly_t L, const fmpq_poly_t A, const fmpq_poly_t B)
    void _fmpq_poly_resultant(fmpz_t rnum, fmpz_t rden, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2)
    void fmpq_poly_resultant(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g)
    void fmpq_poly_resultant_div(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g, const fmpz_t div, slong nbits)
    void _fmpq_poly_derivative(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len)
    void fmpq_poly_derivative(fmpq_poly_t res, const fmpq_poly_t poly)
    void _fmpq_poly_nth_derivative(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, ulong n, slong len)
    void fmpq_poly_nth_derivative(fmpq_poly_t res, const fmpq_poly_t poly, ulong n)
    void _fmpq_poly_integral(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len)
    void fmpq_poly_integral(fmpq_poly_t res, const fmpq_poly_t poly)
    void _fmpq_poly_sqrt_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_sqrt_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_invsqrt_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_invsqrt_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_power_sums(fmpz_struct * res, fmpz_t rden, const fmpz_struct * poly, slong len, slong n)
    void fmpq_poly_power_sums(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void _fmpq_poly_power_sums_to_poly(fmpz_struct * res, const fmpz_struct * poly, const fmpz_t den, slong len)
    void fmpq_poly_power_sums_to_fmpz_poly(fmpz_poly_t res, const fmpq_poly_t Q)
    void fmpq_poly_power_sums_to_poly(fmpq_poly_t res, const fmpq_poly_t Q)
    void _fmpq_poly_log_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_log_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_exp_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * h, const fmpz_t hden, slong hlen, slong n)
    void fmpq_poly_exp_series(fmpq_poly_t res, const fmpq_poly_t h, slong n)
    void _fmpq_poly_exp_expinv_series(fmpz_struct * res1, fmpz_t res1den, fmpz_struct * res2, fmpz_t res2den, const fmpz_struct * h, const fmpz_t hden, slong hlen, slong n)
    void fmpq_poly_exp_expinv_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t h, slong n)
    void _fmpq_poly_atan_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_atan_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_atanh_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_atanh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_asin_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_asin_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_asinh_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_asinh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_tan_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_tan_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_sin_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_sin_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_cos_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_cos_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_sin_cos_series(fmpz_struct * s, fmpz_t sden, fmpz_struct * c, fmpz_t cden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_sin_cos_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t f, slong n)
    void _fmpq_poly_sinh_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_sinh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_cosh_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_cosh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_sinh_cosh_series(fmpz_struct * s, fmpz_t sden, fmpz_struct * c, fmpz_t cden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_sinh_cosh_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t f, slong n)
    void _fmpq_poly_tanh_series(fmpz_struct * g, fmpz_t gden, const fmpz_struct * f, const fmpz_t fden, slong flen, slong n)
    void fmpq_poly_tanh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
    void _fmpq_poly_legendre_p(fmpz_struct * coeffs, fmpz_t den, ulong n)
    void fmpq_poly_legendre_p(fmpq_poly_t poly, ulong n)
    void _fmpq_poly_laguerre_l(fmpz_struct * coeffs, fmpz_t den, ulong n)
    void fmpq_poly_laguerre_l(fmpq_poly_t poly, ulong n)
    void _fmpq_poly_gegenbauer_c(fmpz_struct * coeffs, fmpz_t den, ulong n, const fmpq_t a)
    void fmpq_poly_gegenbauer_c(fmpq_poly_t poly, ulong n, const fmpq_t a)
    void _fmpq_poly_evaluate_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, const fmpz_t a)
    void fmpq_poly_evaluate_fmpz(fmpq_t res, const fmpq_poly_t poly, const fmpz_t a)
    void _fmpq_poly_evaluate_fmpq(fmpz_t rnum, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len, const fmpz_t anum, const fmpz_t aden)
    void fmpq_poly_evaluate_fmpq(fmpq_t res, const fmpq_poly_t poly, const fmpq_t a)
    void _fmpq_poly_interpolate_fmpz_vec(fmpz_struct * poly, fmpz_t den, const fmpz_struct * xs, const fmpz_struct * ys, slong n)
    void fmpq_poly_interpolate_fmpz_vec(fmpq_poly_t poly, const fmpz_struct * xs, const fmpz_struct * ys, slong n)
    void _fmpq_poly_compose(fmpz_struct * res, fmpz_t den, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2)
    void fmpq_poly_compose(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
    void _fmpq_poly_rescale(fmpz_struct * res, fmpz_t denr, const fmpz_struct * poly, const fmpz_t den, slong len, const fmpz_t anum, const fmpz_t aden)
    void fmpq_poly_rescale(fmpq_poly_t res, const fmpq_poly_t poly, const fmpq_t a)
    void _fmpq_poly_compose_series_horner(fmpz_struct * res, fmpz_t den, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, slong n)
    void fmpq_poly_compose_series_horner(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    void _fmpq_poly_compose_series_brent_kung(fmpz_struct * res, fmpz_t den, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, slong n)
    void fmpq_poly_compose_series_brent_kung(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    void _fmpq_poly_compose_series(fmpz_struct * res, fmpz_t den, const fmpz_struct * poly1, const fmpz_t den1, slong len1, const fmpz_struct * poly2, const fmpz_t den2, slong len2, slong n)
    void fmpq_poly_compose_series(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)
    void _fmpq_poly_revert_series_lagrange(fmpz_struct * res, fmpz_t den, const fmpz_struct * poly1, const fmpz_t den1, slong len1, slong n)
    void fmpq_poly_revert_series_lagrange(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void _fmpq_poly_revert_series_lagrange_fast(fmpz_struct * res, fmpz_t den, const fmpz_struct * poly1, const fmpz_t den1, slong len1, slong n)
    void fmpq_poly_revert_series_lagrange_fast(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void _fmpq_poly_revert_series_newton(fmpz_struct * res, fmpz_t den, const fmpz_struct * poly1, const fmpz_t den1, slong len1, slong n)
    void fmpq_poly_revert_series_newton(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void _fmpq_poly_revert_series(fmpz_struct * res, fmpz_t den, const fmpz_struct * poly1, const fmpz_t den1, slong len1, slong n)
    void fmpq_poly_revert_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
    void _fmpq_poly_content(fmpq_t res, const fmpz_struct * poly, const fmpz_t den, slong len)
    void fmpq_poly_content(fmpq_t res, const fmpq_poly_t poly)
    void _fmpq_poly_primitive_part(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len)
    void fmpq_poly_primitive_part(fmpq_poly_t res, const fmpq_poly_t poly)
    int _fmpq_poly_is_monic(const fmpz_struct * poly, const fmpz_t den, slong len)
    int fmpq_poly_is_monic(const fmpq_poly_t poly)
    void _fmpq_poly_make_monic(fmpz_struct * rpoly, fmpz_t rden, const fmpz_struct * poly, const fmpz_t den, slong len)
    void fmpq_poly_make_monic(fmpq_poly_t res, const fmpq_poly_t poly)
    int fmpq_poly_is_squarefree(const fmpq_poly_t poly)
    int _fmpq_poly_print(const fmpz_struct * poly, const fmpz_t den, slong len)
    int fmpq_poly_print(const fmpq_poly_t poly)
    int _fmpq_poly_print_pretty(const fmpz_struct * poly, const fmpz_t den, slong len, const char * x)
    int fmpq_poly_print_pretty(const fmpq_poly_t poly, const char * var)
    # int _fmpq_poly_fprint(FILE * file, const fmpz_struct * poly, const fmpz_t den, slong len)
    # int fmpq_poly_fprint(FILE * file, const fmpq_poly_t poly)
    # int _fmpq_poly_fprint_pretty(FILE * file, const fmpz_struct * poly, const fmpz_t den, slong len, const char * x)
    # int fmpq_poly_fprint_pretty(FILE * file, const fmpq_poly_t poly, const char * var)
    int fmpq_poly_read(fmpq_poly_t poly)
    # int fmpq_poly_fread(FILE * file, fmpq_poly_t poly)
