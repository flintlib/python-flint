from flint.flintlib.fq_default cimport fq_default_ctx_t, fq_default_t
from flint.flintlib.nmod_poly cimport nmod_poly_t
from flint.flintlib.flint cimport slong, ulong, flint_rand_t
from flint.flintlib.fmpz_poly cimport fmpz_poly_t
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.fmpz_mod_poly cimport fmpz_mod_poly_t
from flint.flintlib.fq_poly cimport fq_poly_t
from flint.flintlib.fq_nmod_poly cimport fq_nmod_poly_t
from flint.flintlib.fq_zech_poly cimport fq_zech_poly_t


cdef extern from "flint/fq_default_poly.h":
    # Type definitions **********************************************/
    ctypedef union fq_default_poly_struct:
        fq_poly_t fq
        fq_nmod_poly_t fq_nmod
        fq_zech_poly_t fq_zech
        nmod_poly_t nmod
        fmpz_mod_poly_t fmpz_mod
    ctypedef fq_default_poly_struct fq_default_poly_t[1]

    # Parsed from here **********************************************/
    void fq_default_poly_init(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    void fq_default_poly_init2(fq_default_poly_t poly, slong alloc, const fq_default_ctx_t ctx)
    void fq_default_poly_realloc(fq_default_poly_t poly, slong alloc, const fq_default_ctx_t ctx)
    void fq_default_poly_fit_length(fq_default_poly_t poly, slong len, const fq_default_ctx_t ctx)
    void fq_default_poly_clear(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    void _fq_default_poly_set_length(fq_default_poly_t poly, slong len, const fq_default_ctx_t ctx)
    void fq_default_poly_truncate(fq_default_poly_t poly, slong newlen, const fq_default_ctx_t ctx)
    void fq_default_poly_set_trunc(fq_default_poly_t poly1, fq_default_poly_t poly2, slong newlen, const fq_default_ctx_t ctx)
    void fq_default_poly_reverse(fq_default_poly_t output, const fq_default_poly_t input, slong m, const fq_default_ctx_t ctx)
    slong fq_default_poly_degree(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    slong fq_default_poly_length(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    void fq_default_poly_randtest(fq_default_poly_t f, flint_rand_t state, slong len, const fq_default_ctx_t ctx)
    void fq_default_poly_randtest_not_zero(fq_default_poly_t f, flint_rand_t state, slong len, const fq_default_ctx_t ctx)
    void fq_default_poly_randtest_monic(fq_default_poly_t f, flint_rand_t state, slong len, const fq_default_ctx_t ctx)
    void fq_default_poly_randtest_irreducible(fq_default_poly_t f, flint_rand_t state, slong len, const fq_default_ctx_t ctx)
    void fq_default_poly_set(fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
    void fq_default_poly_set_fq_default(fq_default_poly_t poly, const fq_default_t c, const fq_default_ctx_t ctx)
    void fq_default_poly_swap(fq_default_poly_t op1, fq_default_poly_t op2, const fq_default_ctx_t ctx)
    void fq_default_poly_zero(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    void fq_default_poly_one(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    void fq_default_poly_gen(fq_default_poly_t poly, const fq_default_ctx_t ctx)
    void fq_default_poly_make_monic(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)
    void fq_default_poly_set_nmod_poly(fq_default_poly_t rop, const nmod_poly_t op, const fq_default_ctx_t ctx)
    void fq_default_poly_set_fmpz_mod_poly(fq_default_poly_t rop, const fmpz_mod_poly_t op, const fq_default_ctx_t ctx)
    void fq_default_poly_set_fmpz_poly(fq_default_poly_t rop, const fmpz_poly_t op, const fq_default_ctx_t ctx)
    void fq_default_poly_get_coeff(fq_default_t x, const fq_default_poly_t poly, slong n, const fq_default_ctx_t ctx)
    void fq_default_poly_set_coeff(fq_default_poly_t poly, slong n, const fq_default_t x, const fq_default_ctx_t ctx)
    void fq_default_poly_set_coeff_fmpz(fq_default_poly_t poly, slong n, const fmpz_t x, const fq_default_ctx_t ctx)
    int fq_default_poly_equal(const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
    int fq_default_poly_equal_trunc(const fq_default_poly_t poly1, const fq_default_poly_t poly2, slong n, const fq_default_ctx_t ctx)
    int fq_default_poly_is_zero(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    int fq_default_poly_is_one(const fq_default_poly_t op, const fq_default_ctx_t ctx)
    int fq_default_poly_is_gen(const fq_default_poly_t op, const fq_default_ctx_t ctx)
    int fq_default_poly_is_unit(const fq_default_poly_t op, const fq_default_ctx_t ctx)
    int fq_default_poly_equal_fq_default(const fq_default_poly_t poly, const fq_default_t c, const fq_default_ctx_t ctx)
    void fq_default_poly_add(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
    void fq_default_poly_add_si(fq_default_poly_t res, const fq_default_poly_t poly1, slong c, const fq_default_ctx_t ctx)
    void fq_default_poly_add_series(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, slong n, const fq_default_ctx_t ctx)
    void fq_default_poly_sub(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
    void fq_default_poly_sub_series(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, slong n, const fq_default_ctx_t ctx)
    void fq_default_poly_neg(fq_default_poly_t res, const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    void fq_default_poly_scalar_mul_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
    void fq_default_poly_scalar_addmul_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
    void fq_default_poly_scalar_submul_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
    void fq_default_poly_scalar_div_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
    void fq_default_poly_mul(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)
    void fq_default_poly_mullow(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, slong n, const fq_default_ctx_t ctx)
    void fq_default_poly_mulhigh(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, slong start, const fq_default_ctx_t ctx)
    void fq_default_poly_mulmod(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_poly_t f, const fq_default_ctx_t ctx)
    void fq_default_poly_sqr(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)
    void fq_default_poly_pow(fq_default_poly_t rop, const fq_default_poly_t op, ulong e, const fq_default_ctx_t ctx)
    void fq_default_poly_powmod_ui_binexp(fq_default_poly_t res, const fq_default_poly_t poly, ulong e, const fq_default_poly_t f, const fq_default_ctx_t ctx)
    void fq_default_poly_powmod_fmpz_binexp(fq_default_poly_t res, const fq_default_poly_t poly, const fmpz_t e, const fq_default_poly_t f, const fq_default_ctx_t ctx)
    void fq_default_poly_pow_trunc(fq_default_poly_t res, const fq_default_poly_t poly, ulong e, slong trunc, const fq_default_ctx_t ctx)
    void fq_default_poly_shift_left(fq_default_poly_t rop, const fq_default_poly_t op, slong n, const fq_default_ctx_t ctx)
    void fq_default_poly_shift_right(fq_default_poly_t rop, const fq_default_poly_t op, slong n, const fq_default_ctx_t ctx)
    slong fq_default_poly_hamming_weight(const fq_default_poly_t op, const fq_default_ctx_t ctx)
    void fq_default_poly_divrem(fq_default_poly_t Q, fq_default_poly_t R, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)
    void fq_default_poly_rem(fq_default_poly_t R, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)
    void fq_default_poly_inv_series(fq_default_poly_t Qinv, const fq_default_poly_t Q, slong n, const fq_default_ctx_t ctx)
    void fq_default_poly_div_series(fq_default_poly_t Q, const fq_default_poly_t A, const fq_default_poly_t B, slong n, const fq_default_ctx_t ctx)
    void fq_default_poly_gcd(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)
    void fq_default_poly_xgcd(fq_default_poly_t G, fq_default_poly_t S, fq_default_poly_t T, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)
    int fq_default_poly_divides(fq_default_poly_t Q, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)
    void fq_default_poly_derivative(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)
    void fq_default_poly_invsqrt_series(fq_default_poly_t g, const fq_default_poly_t h, slong n, fq_default_ctx_t ctx)
    void fq_default_poly_sqrt_series(fq_default_poly_t g, const fq_default_poly_t h, slong n, fq_default_ctx_t ctx)
    int fq_default_poly_sqrt(fq_default_poly_t s, const fq_default_poly_t p, fq_default_ctx_t mod)
    void fq_default_poly_evaluate_fq_default(fq_default_t rop, const fq_default_poly_t f, const fq_default_t a, const fq_default_ctx_t ctx)
    void fq_default_poly_compose(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)
    void fq_default_poly_compose_mod(fq_default_poly_t res, const fq_default_poly_t f, const fq_default_poly_t g, const fq_default_poly_t h, const fq_default_ctx_t ctx)
    # int fq_default_poly_fprint_pretty(FILE * file, const fq_default_poly_t poly, const char * x, const fq_default_ctx_t ctx)
    int fq_default_poly_print_pretty(const fq_default_poly_t poly, const char * x, const fq_default_ctx_t ctx)
    # int fq_default_poly_fprint(FILE * file, const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    int fq_default_poly_print(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    char * fq_default_poly_get_str(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
    char * fq_default_poly_get_str_pretty(const fq_default_poly_t poly, const char * x, const fq_default_ctx_t ctx)
    void fq_default_poly_inflate(fq_default_poly_t result, const fq_default_poly_t input, ulong inflation, const fq_default_ctx_t ctx)
    void fq_default_poly_deflate(fq_default_poly_t result, const fq_default_poly_t input, ulong deflation, const fq_default_ctx_t ctx)
    ulong fq_default_poly_deflation(const fq_default_poly_t input, const fq_default_ctx_t ctx)