from flint.flintlib.types.acb cimport acb_srcptr, acb_t
from flint.flintlib.types.arb cimport arb_poly_t, arb_ptr, arb_srcptr, arb_t, mag_t
from flint.flintlib.types.arf cimport arf_t
from flint.flintlib.types.flint cimport flint_rand_t, slong, ulong
from flint.flintlib.types.fmpq cimport fmpq_poly_t
from flint.flintlib.types.fmpz cimport fmpz_poly_t

# unknown type FILE

# .. macro:: arb_poly_get_coeff_ptr(poly, n)

cdef extern from "flint/arb_poly.h":
    void arb_poly_init(arb_poly_t poly)
    void arb_poly_clear(arb_poly_t poly)
    void arb_poly_fit_length(arb_poly_t poly, slong len)
    void _arb_poly_set_length(arb_poly_t poly, slong len)
    void _arb_poly_normalise(arb_poly_t poly)
    slong arb_poly_allocated_bytes(const arb_poly_t x)
    slong arb_poly_length(const arb_poly_t poly)
    slong arb_poly_degree(const arb_poly_t poly)
    int arb_poly_is_zero(const arb_poly_t poly)
    int arb_poly_is_one(const arb_poly_t poly)
    int arb_poly_is_x(const arb_poly_t poly)
    void arb_poly_zero(arb_poly_t poly)
    void arb_poly_one(arb_poly_t poly)
    void arb_poly_set(arb_poly_t dest, const arb_poly_t src)
    void arb_poly_set_round(arb_poly_t dest, const arb_poly_t src, slong prec)
    void arb_poly_set_trunc(arb_poly_t dest, const arb_poly_t src, slong n)
    void arb_poly_set_trunc_round(arb_poly_t dest, const arb_poly_t src, slong n, slong prec)
    void arb_poly_set_coeff_si(arb_poly_t poly, slong n, slong c)
    void arb_poly_set_coeff_arb(arb_poly_t poly, slong n, const arb_t c)
    void arb_poly_get_coeff_arb(arb_t v, const arb_poly_t poly, slong n)
    void _arb_poly_shift_right(arb_ptr res, arb_srcptr poly, slong len, slong n)
    void arb_poly_shift_right(arb_poly_t res, const arb_poly_t poly, slong n)
    void _arb_poly_shift_left(arb_ptr res, arb_srcptr poly, slong len, slong n)
    void arb_poly_shift_left(arb_poly_t res, const arb_poly_t poly, slong n)
    void arb_poly_truncate(arb_poly_t poly, slong n)
    slong arb_poly_valuation(const arb_poly_t poly)
    void arb_poly_set_fmpz_poly(arb_poly_t poly, const fmpz_poly_t src, slong prec)
    void arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src, slong prec)
    void arb_poly_set_si(arb_poly_t poly, slong src)
    void arb_poly_printd(const arb_poly_t poly, slong digits)
    # void arb_poly_fprintd(FILE * file, const arb_poly_t poly, slong digits)
    void arb_poly_randtest(arb_poly_t poly, flint_rand_t state, slong len, slong prec, slong mag_bits)
    int arb_poly_contains(const arb_poly_t poly1, const arb_poly_t poly2)
    int arb_poly_contains_fmpz_poly(const arb_poly_t poly1, const fmpz_poly_t poly2)
    int arb_poly_contains_fmpq_poly(const arb_poly_t poly1, const fmpq_poly_t poly2)
    int arb_poly_equal(const arb_poly_t A, const arb_poly_t B)
    int _arb_poly_overlaps(arb_srcptr poly1, slong len1, arb_srcptr poly2, slong len2)
    int arb_poly_overlaps(const arb_poly_t poly1, const arb_poly_t poly2)
    int arb_poly_get_unique_fmpz_poly(fmpz_poly_t z, const arb_poly_t x)
    void _arb_poly_majorant(arb_ptr res, arb_srcptr poly, slong len, slong prec)
    void arb_poly_majorant(arb_poly_t res, const arb_poly_t poly, slong prec)
    void _arb_poly_add(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
    void arb_poly_add(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong prec)
    void arb_poly_add_si(arb_poly_t C, const arb_poly_t A, slong B, slong prec)
    void _arb_poly_sub(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
    void arb_poly_sub(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong prec)
    void arb_poly_add_series(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong len, slong prec)
    void arb_poly_sub_series(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong len, slong prec)
    void arb_poly_neg(arb_poly_t C, const arb_poly_t A)
    void arb_poly_scalar_mul_2exp_si(arb_poly_t C, const arb_poly_t A, slong c)
    void arb_poly_scalar_mul(arb_poly_t C, const arb_poly_t A, const arb_t c, slong prec)
    void arb_poly_scalar_div(arb_poly_t C, const arb_poly_t A, const arb_t c, slong prec)
    void _arb_poly_mullow_classical(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong n, slong prec)
    void _arb_poly_mullow_block(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong n, slong prec)
    void _arb_poly_mullow(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong n, slong prec)
    void arb_poly_mullow_classical(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong n, slong prec)
    void arb_poly_mullow_block(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong n, slong prec)
    void arb_poly_mullow(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong n, slong prec)
    void _arb_poly_mul(arb_ptr C, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
    void arb_poly_mul(arb_poly_t C, const arb_poly_t A, const arb_poly_t B, slong prec)
    void _arb_poly_inv_series(arb_ptr Q, arb_srcptr A, slong Alen, slong len, slong prec)
    void arb_poly_inv_series(arb_poly_t Q, const arb_poly_t A, slong n, slong prec)
    void _arb_poly_div_series(arb_ptr Q, arb_srcptr A, slong Alen, arb_srcptr B, slong Blen, slong n, slong prec)
    void arb_poly_div_series(arb_poly_t Q, const arb_poly_t A, const arb_poly_t B, slong n, slong prec)
    void _arb_poly_div(arb_ptr Q, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
    void _arb_poly_rem(arb_ptr R, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
    void _arb_poly_divrem(arb_ptr Q, arb_ptr R, arb_srcptr A, slong lenA, arb_srcptr B, slong lenB, slong prec)
    int arb_poly_divrem(arb_poly_t Q, arb_poly_t R, const arb_poly_t A, const arb_poly_t B, slong prec)
    void _arb_poly_div_root(arb_ptr Q, arb_t R, arb_srcptr A, slong len, const arb_t c, slong prec)
    void _arb_poly_taylor_shift(arb_ptr g, const arb_t c, slong n, slong prec)
    void arb_poly_taylor_shift(arb_poly_t g, const arb_poly_t f, const arb_t c, slong prec)
    void _arb_poly_compose(arb_ptr res, arb_srcptr poly1, slong len1, arb_srcptr poly2, slong len2, slong prec)
    void arb_poly_compose(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, slong prec)
    void _arb_poly_compose_series(arb_ptr res, arb_srcptr poly1, slong len1, arb_srcptr poly2, slong len2, slong n, slong prec)
    void arb_poly_compose_series(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, slong n, slong prec)
    void _arb_poly_revert_series(arb_ptr h, arb_srcptr f, slong flen, slong n, slong prec)
    void arb_poly_revert_series(arb_poly_t h, const arb_poly_t f, slong n, slong prec)
    void _arb_poly_evaluate_horner(arb_t y, arb_srcptr f, slong len, const arb_t x, slong prec)
    void arb_poly_evaluate_horner(arb_t y, const arb_poly_t f, const arb_t x, slong prec)
    void _arb_poly_evaluate_rectangular(arb_t y, arb_srcptr f, slong len, const arb_t x, slong prec)
    void arb_poly_evaluate_rectangular(arb_t y, const arb_poly_t f, const arb_t x, slong prec)
    void _arb_poly_evaluate(arb_t y, arb_srcptr f, slong len, const arb_t x, slong prec)
    void arb_poly_evaluate(arb_t y, const arb_poly_t f, const arb_t x, slong prec)
    void _arb_poly_evaluate_acb_horner(acb_t y, arb_srcptr f, slong len, const acb_t x, slong prec)
    void arb_poly_evaluate_acb_horner(acb_t y, const arb_poly_t f, const acb_t x, slong prec)
    void _arb_poly_evaluate_acb_rectangular(acb_t y, arb_srcptr f, slong len, const acb_t x, slong prec)
    void arb_poly_evaluate_acb_rectangular(acb_t y, const arb_poly_t f, const acb_t x, slong prec)
    void _arb_poly_evaluate_acb(acb_t y, arb_srcptr f, slong len, const acb_t x, slong prec)
    void arb_poly_evaluate_acb(acb_t y, const arb_poly_t f, const acb_t x, slong prec)
    void _arb_poly_evaluate2_horner(arb_t y, arb_t z, arb_srcptr f, slong len, const arb_t x, slong prec)
    void arb_poly_evaluate2_horner(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, slong prec)
    void _arb_poly_evaluate2_rectangular(arb_t y, arb_t z, arb_srcptr f, slong len, const arb_t x, slong prec)
    void arb_poly_evaluate2_rectangular(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, slong prec)
    void _arb_poly_evaluate2(arb_t y, arb_t z, arb_srcptr f, slong len, const arb_t x, slong prec)
    void arb_poly_evaluate2(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, slong prec)
    void _arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec)
    void arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, slong prec)
    void _arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec)
    void arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, slong prec)
    void _arb_poly_evaluate2_acb(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec)
    void arb_poly_evaluate2_acb(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, slong prec)
    void _arb_poly_product_roots(arb_ptr poly, arb_srcptr xs, slong n, slong prec)
    void arb_poly_product_roots(arb_poly_t poly, arb_srcptr xs, slong n, slong prec)
    void _arb_poly_product_roots_complex(arb_ptr poly, arb_srcptr r, slong rn, acb_srcptr c, slong cn, slong prec)
    void arb_poly_product_roots_complex(arb_poly_t poly, arb_srcptr r, slong rn, acb_srcptr c, slong cn, slong prec)
    arb_ptr * _arb_poly_tree_alloc(slong len)
    void _arb_poly_tree_free(arb_ptr * tree, slong len)
    void _arb_poly_tree_build(arb_ptr * tree, arb_srcptr roots, slong len, slong prec)
    void _arb_poly_evaluate_vec_iter(arb_ptr ys, arb_srcptr poly, slong plen, arb_srcptr xs, slong n, slong prec)
    void arb_poly_evaluate_vec_iter(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, slong n, slong prec)
    void _arb_poly_evaluate_vec_fast_precomp(arb_ptr vs, arb_srcptr poly, slong plen, arb_ptr * tree, slong len, slong prec)
    void _arb_poly_evaluate_vec_fast(arb_ptr ys, arb_srcptr poly, slong plen, arb_srcptr xs, slong n, slong prec)
    void arb_poly_evaluate_vec_fast(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, slong n, slong prec)
    void _arb_poly_interpolate_newton(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
    void arb_poly_interpolate_newton(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
    void _arb_poly_interpolate_barycentric(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
    void arb_poly_interpolate_barycentric(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
    void _arb_poly_interpolation_weights(arb_ptr w, arb_ptr * tree, slong len, slong prec)
    void _arb_poly_interpolate_fast_precomp(arb_ptr poly, arb_srcptr ys, arb_ptr * tree, arb_srcptr weights, slong len, slong prec)
    void _arb_poly_interpolate_fast(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, slong len, slong prec)
    void arb_poly_interpolate_fast(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, slong n, slong prec)
    void _arb_poly_derivative(arb_ptr res, arb_srcptr poly, slong len, slong prec)
    void arb_poly_derivative(arb_poly_t res, const arb_poly_t poly, slong prec)
    void _arb_poly_nth_derivative(arb_ptr res, arb_srcptr poly, ulong n, slong len, slong prec)
    void arb_poly_nth_derivative(arb_poly_t res, const arb_poly_t poly, ulong n, slong prec)
    void _arb_poly_integral(arb_ptr res, arb_srcptr poly, slong len, slong prec)
    void arb_poly_integral(arb_poly_t res, const arb_poly_t poly, slong prec)
    void _arb_poly_borel_transform(arb_ptr res, arb_srcptr poly, slong len, slong prec)
    void arb_poly_borel_transform(arb_poly_t res, const arb_poly_t poly, slong prec)
    void _arb_poly_inv_borel_transform(arb_ptr res, arb_srcptr poly, slong len, slong prec)
    void arb_poly_inv_borel_transform(arb_poly_t res, const arb_poly_t poly, slong prec)
    void _arb_poly_binomial_transform_basecase(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec)
    void arb_poly_binomial_transform_basecase(arb_poly_t b, const arb_poly_t a, slong len, slong prec)
    void _arb_poly_binomial_transform_convolution(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec)
    void arb_poly_binomial_transform_convolution(arb_poly_t b, const arb_poly_t a, slong len, slong prec)
    void _arb_poly_binomial_transform(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec)
    void arb_poly_binomial_transform(arb_poly_t b, const arb_poly_t a, slong len, slong prec)
    void _arb_poly_graeffe_transform(arb_ptr b, arb_srcptr a, slong len, slong prec)
    void arb_poly_graeffe_transform(arb_poly_t b, const arb_poly_t a, slong prec)
    void _arb_poly_pow_ui_trunc_binexp(arb_ptr res, arb_srcptr f, slong flen, ulong exp, slong len, slong prec)
    void arb_poly_pow_ui_trunc_binexp(arb_poly_t res, const arb_poly_t poly, ulong exp, slong len, slong prec)
    void _arb_poly_pow_ui(arb_ptr res, arb_srcptr f, slong flen, ulong exp, slong prec)
    void arb_poly_pow_ui(arb_poly_t res, const arb_poly_t poly, ulong exp, slong prec)
    void _arb_poly_pow_series(arb_ptr h, arb_srcptr f, slong flen, arb_srcptr g, slong glen, slong len, slong prec)
    void arb_poly_pow_series(arb_poly_t h, const arb_poly_t f, const arb_poly_t g, slong len, slong prec)
    void _arb_poly_pow_arb_series(arb_ptr h, arb_srcptr f, slong flen, const arb_t g, slong len, slong prec)
    void arb_poly_pow_arb_series(arb_poly_t h, const arb_poly_t f, const arb_t g, slong len, slong prec)
    void _arb_poly_sqrt_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sqrt_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_rsqrt_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_rsqrt_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_log_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
    void arb_poly_log_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
    void _arb_poly_log1p_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
    void arb_poly_log1p_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
    void _arb_poly_atan_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
    void arb_poly_atan_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
    void _arb_poly_asin_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
    void arb_poly_asin_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
    void _arb_poly_acos_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec)
    void arb_poly_acos_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec)
    void _arb_poly_exp_series_basecase(arb_ptr f, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_exp_series_basecase(arb_poly_t f, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_exp_series(arb_ptr f, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_exp_series(arb_poly_t f, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sin_cos_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sin_cos_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sin_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sin_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_cos_series(arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_cos_series(arb_poly_t c, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_tan_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec)
    void arb_poly_tan_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sin_cos_pi_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sin_cos_pi_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sin_pi_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sin_pi_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_cos_pi_series(arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_cos_pi_series(arb_poly_t c, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_cot_pi_series(arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_cot_pi_series(arb_poly_t c, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sinh_cosh_series_basecase(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sinh_cosh_series_basecase(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sinh_cosh_series_exponential(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sinh_cosh_series_exponential(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sinh_cosh_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sinh_cosh_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sinh_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sinh_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_cosh_series(arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_cosh_series(arb_poly_t c, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sinc_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sinc_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_sinc_pi_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_sinc_pi_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_lambertw_series(arb_ptr res, arb_srcptr z, slong zlen, int flags, slong len, slong prec)
    void arb_poly_lambertw_series(arb_poly_t res, const arb_poly_t z, int flags, slong len, slong prec)
    void _arb_poly_gamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_gamma_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_rgamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_rgamma_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_lgamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_lgamma_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_digamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_digamma_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_rising_ui_series(arb_ptr res, arb_srcptr f, slong flen, ulong r, slong trunc, slong prec)
    void arb_poly_rising_ui_series(arb_poly_t res, const arb_poly_t f, ulong r, slong trunc, slong prec)
    void arb_poly_zeta_series(arb_poly_t res, const arb_poly_t s, const arb_t a, int deflate, slong n, slong prec)
    void _arb_poly_riemann_siegel_theta_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_riemann_siegel_theta_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_riemann_siegel_z_series(arb_ptr res, arb_srcptr h, slong hlen, slong n, slong prec)
    void arb_poly_riemann_siegel_z_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec)
    void _arb_poly_root_bound_fujiwara(mag_t bound, arb_srcptr poly, slong len)
    void arb_poly_root_bound_fujiwara(mag_t bound, arb_poly_t poly)
    void _arb_poly_newton_convergence_factor(arf_t convergence_factor, arb_srcptr poly, slong len, const arb_t convergence_interval, slong prec)
    int _arb_poly_newton_step(arb_t xnew, arb_srcptr poly, slong len, const arb_t x, const arb_t convergence_interval, const arf_t convergence_factor, slong prec)
    void _arb_poly_newton_refine_root(arb_t r, arb_srcptr poly, slong len, const arb_t start, const arb_t convergence_interval, const arf_t convergence_factor, slong eval_extra_prec, slong prec)
    void _arb_poly_swinnerton_dyer_ui(arb_ptr poly, ulong n, slong trunc, slong prec)
    void arb_poly_swinnerton_dyer_ui(arb_poly_t poly, ulong n, slong prec)
