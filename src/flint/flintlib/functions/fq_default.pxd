from flint.flintlib.types.flint cimport flint_rand_t, fmpz_t, slong, ulong
from flint.flintlib.types.fmpz cimport fmpz_poly_t
from flint.flintlib.types.fmpz_mod cimport fmpz_mod_ctx_t
from flint.flintlib.types.fmpz_mod_poly cimport fmpz_mod_poly_t
from flint.flintlib.types.fq_default cimport fq_default_ctx_t, fq_default_t
from flint.flintlib.types.nmod cimport nmod_poly_t

# unknown type FILE


cdef extern from "flint/fq_default.h":
    void fq_default_ctx_init_type(fq_default_ctx_t ctx, const fmpz_t p, slong d, const char * var, int type)
    void fq_default_ctx_init(fq_default_ctx_t ctx, const fmpz_t p, slong d, const char * var)
    void fq_default_ctx_init_modulus_nmod_type(fq_default_ctx_t ctx, const nmod_poly_t modulus, const char * var, int type)
    void fq_default_ctx_init_modulus_nmod(fq_default_ctx_t ctx, const nmod_poly_t modulus, const char * var)
    void fq_default_ctx_init_modulus_type(fq_default_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var, int type)
    void fq_default_ctx_init_modulus(fq_default_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var)
    void fq_default_ctx_clear(fq_default_ctx_t ctx)
    int fq_default_ctx_type(const fq_default_ctx_t ctx)
    void * fq_default_ctx_inner(const fq_default_ctx_t ctx)
    slong fq_default_ctx_degree(const fq_default_ctx_t ctx)
    void fq_default_ctx_prime(fmpz_t prime, const fq_default_ctx_t ctx)
    void fq_default_ctx_order(fmpz_t f, const fq_default_ctx_t ctx)
    void fq_default_ctx_modulus(fmpz_mod_poly_t p, const fq_default_ctx_t ctx)
    # int fq_default_ctx_fprint(FILE * file, const fq_default_ctx_t ctx)
    void fq_default_ctx_print(const fq_default_ctx_t ctx)
    void fq_default_ctx_init_randtest(fq_default_ctx_t ctx, flint_rand_t state)
    void fq_default_get_coeff_fmpz(fmpz_t c, fq_default_t op, slong n, const fq_default_ctx_t ctx)
    void fq_default_init(fq_default_t rop, const fq_default_ctx_t ctx)
    void fq_default_init2(fq_default_t rop, const fq_default_ctx_t ctx)
    void fq_default_clear(fq_default_t rop, const fq_default_ctx_t ctx)
    int fq_default_is_invertible(const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_add(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)
    void fq_default_sub(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)
    void fq_default_sub_one(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)
    void fq_default_neg(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_mul(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)
    void fq_default_mul_fmpz(fq_default_t rop, const fq_default_t op, const fmpz_t x, const fq_default_ctx_t ctx)
    void fq_default_mul_si(fq_default_t rop, const fq_default_t op, slong x, const fq_default_ctx_t ctx)
    void fq_default_mul_ui(fq_default_t rop, const fq_default_t op, ulong x, const fq_default_ctx_t ctx)
    void fq_default_sqr(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_div(fq_default_t rop, fq_default_t op1, fq_default_t op2, const fq_default_ctx_t ctx)
    void fq_default_inv(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_pow(fq_default_t rop, const fq_default_t op, const fmpz_t e, const fq_default_ctx_t ctx)
    void fq_default_pow_ui(fq_default_t rop, const fq_default_t op, const ulong e, const fq_default_ctx_t ctx)
    int fq_default_sqrt(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)
    void fq_default_pth_root(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)
    int fq_default_is_square(const fq_default_t op, const fq_default_ctx_t ctx)
    # int fq_default_fprint_pretty(FILE * file, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_print_pretty(const fq_default_t op, const fq_default_ctx_t ctx)
    # int fq_default_fprint(FILE * file, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_print(const fq_default_t op, const fq_default_ctx_t ctx)
    char * fq_default_get_str(const fq_default_t op, const fq_default_ctx_t ctx)
    char * fq_default_get_str_pretty(const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_randtest(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)
    void fq_default_randtest_not_zero(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)
    void fq_default_rand(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)
    void fq_default_rand_not_zero(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)
    void fq_default_set(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_set_si(fq_default_t rop, const slong x, const fq_default_ctx_t ctx)
    void fq_default_set_ui(fq_default_t rop, const ulong x, const fq_default_ctx_t ctx)
    void fq_default_set_fmpz(fq_default_t rop, const fmpz_t x, const fq_default_ctx_t ctx)
    void fq_default_swap(fq_default_t op1, fq_default_t op2, const fq_default_ctx_t ctx)
    void fq_default_zero(fq_default_t rop, const fq_default_ctx_t ctx)
    void fq_default_one(fq_default_t rop, const fq_default_ctx_t ctx)
    void fq_default_gen(fq_default_t rop, const fq_default_ctx_t ctx)
    int fq_default_get_fmpz(fmpz_t rop, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_get_nmod_poly(nmod_poly_t poly, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_set_nmod_poly(fq_default_t op, const nmod_poly_t poly, const fq_default_ctx_t ctx)
    void fq_default_get_fmpz_mod_poly(fmpz_mod_poly_t poly, const fq_default_t op,  const fq_default_ctx_t ctx)
    void fq_default_set_fmpz_mod_poly(fq_default_t op, const fmpz_mod_poly_t poly, const fq_default_ctx_t ctx)
    void fq_default_get_fmpz_poly(fmpz_poly_t a, const fq_default_t b, const fq_default_ctx_t ctx)
    void fq_default_set_fmpz_poly(fq_default_t a, const fmpz_poly_t b, const fq_default_ctx_t ctx)
    int fq_default_is_zero(const fq_default_t op, const fq_default_ctx_t ctx)
    int fq_default_is_one(const fq_default_t op, const fq_default_ctx_t ctx)
    int fq_default_equal(const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)
    void fq_default_trace(fmpz_t rop, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_norm(fmpz_t rop, const fq_default_t op, const fq_default_ctx_t ctx)
    void fq_default_frobenius(fq_default_t rop, const fq_default_t op, slong e, const fq_default_ctx_t ctx)
