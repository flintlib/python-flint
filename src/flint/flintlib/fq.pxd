from flint.flintlib.flint cimport flint_bitcnt_t, fmpz_struct, slong, flint_rand_t, ulong
from flint.flintlib.fmpz cimport fmpz_t, fmpz_struct
from flint.flintlib.fmpz_mod cimport fmpz_mod_ctx_t
from flint.flintlib.fmpz_mod_mat cimport fmpz_mod_mat_t
from flint.flintlib.fmpz_poly cimport fmpz_poly_t, fmpz_poly_struct
from flint.flintlib.fmpz_mod_poly cimport fmpz_mod_poly_t, fmpz_mod_poly_struct

cdef extern from "flint/fq.h":

    ctypedef fmpz_poly_t fq_t
    ctypedef fmpz_poly_struct fq_struct

    ctypedef struct fq_ctx_struct:
        fmpz_mod_ctx_t ctxp

        int sparse_modulus
        int is_conway # whether field was initialized with the Flint Conway tables  (assures primitivity)

        fmpz_struct * a
        slong * j
        slong len

        fmpz_mod_poly_t modulus
        fmpz_mod_poly_t inv

        char * var

    ctypedef fq_ctx_struct fq_ctx_t[1]

    void fq_ctx_init(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    int _fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    void fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    void fq_ctx_init_modulus(fq_ctx_t ctx, const fmpz_mod_poly_t modulus, const fmpz_mod_ctx_t ctxp, const char *var)
    void fq_ctx_clear(fq_ctx_t ctx)
    const fmpz_mod_poly_struct* fq_ctx_modulus(const fq_ctx_t ctx)
    long fq_ctx_degree(const fq_ctx_t ctx)
    fmpz_struct * fq_ctx_prime(const fq_ctx_t ctx)
    void fq_ctx_order(fmpz_t f, const fq_ctx_t ctx)
    # int fq_ctx_fprint(FILE * file, const fq_ctx_t ctx)
    void fq_ctx_print(const fq_ctx_t ctx)
    void fq_ctx_randtest(fq_ctx_t ctx)
    void fq_ctx_randtest_reducible(fq_ctx_t ctx)
    void fq_init(fq_t rop, const fq_ctx_t ctx)
    void fq_init2(fq_t rop, const fq_ctx_t ctx)
    void fq_clear(fq_t rop, const fq_ctx_t ctx)
    void _fq_sparse_reduce(fmpz_struct *R, slong lenR, const fq_ctx_t ctx)
    void _fq_dense_reduce(fmpz_struct *R, slong lenR, const fq_ctx_t ctx)
    void _fq_reduce(fmpz_struct *r, slong lenR, const fq_ctx_t ctx)
    void fq_reduce(fq_t rop, const fq_ctx_t ctx)
    void fq_add(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    void fq_sub(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    void fq_sub_one(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
    void fq_neg(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    void fq_mul(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    void fq_mul_fmpz(fq_t rop, const fq_t op, const fmpz_t x, const fq_ctx_t ctx)
    void fq_mul_si(fq_t rop, const fq_t op, slong x, const fq_ctx_t ctx)
    void fq_mul_ui(fq_t rop, const fq_t op, ulong x, const fq_ctx_t ctx)
    void fq_sqr(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    void fq_div(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    void _fq_inv(fmpz_struct *rop, const fmpz_struct *op, slong len, const fq_ctx_t ctx)
    void fq_inv(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    void fq_gcdinv(fq_t f, fq_t inv, const fq_t op, const fq_ctx_t ctx)
    void _fq_pow(fmpz_struct *rop, const fmpz_struct *op, slong len, const fmpz_t e, const fq_ctx_t ctx)
    void fq_pow(fq_t rop, const fq_t op, const fmpz_t e, const fq_ctx_t ctx)
    void fq_pow_ui(fq_t rop, const fq_t op, const ulong e, const fq_ctx_t ctx)
    int fq_sqrt(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
    void fq_pth_root(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
    int fq_is_square(const fq_t op, const fq_ctx_t ctx)
    # int fq_fprint_pretty(FILE *file, const fq_t op, const fq_ctx_t ctx)
    int fq_print_pretty(const fq_t op, const fq_ctx_t ctx)
    # void fq_fprint(FILE * file, const fq_t op, const fq_ctx_t ctx)
    void fq_print(const fq_t op, const fq_ctx_t ctx)
    char * fq_get_str(const fq_t op, const fq_ctx_t ctx)
    char * fq_get_str_pretty(const fq_t op, const fq_ctx_t ctx)
    void fq_randtest(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    void fq_randtest_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    void fq_randtest_dense(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    void fq_rand(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    void fq_rand_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
    void fq_set(fq_t rop, const fq_t op, const fq_ctx_t ctx)
    void fq_set_si(fq_t rop, const slong x, const fq_ctx_t ctx)
    void fq_set_ui(fq_t rop, const ulong x, const fq_ctx_t ctx)
    void fq_set_fmpz(fq_t rop, const fmpz_t x, const fq_ctx_t ctx)
    void fq_swap(fq_t op1, fq_t op2, const fq_ctx_t ctx)
    void fq_zero(fq_t rop, const fq_ctx_t ctx)
    void fq_one(fq_t rop, const fq_ctx_t ctx)
    void fq_gen(fq_t rop, const fq_ctx_t ctx)
    int fq_get_fmpz(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
    void fq_get_fmpz_poly(fmpz_poly_t a, const fq_t b, const fq_ctx_t ctx)
    void fq_get_fmpz_mod_poly(fmpz_mod_poly_t a, const fq_t b, const fq_ctx_t ctx)
    void fq_set_fmpz_poly(fq_t a, const fmpz_poly_t b, const fq_ctx_t ctx)
    void fq_set_fmpz_mod_poly(fq_t a, const fmpz_mod_poly_t b, const fq_ctx_t ctx)
    void fq_get_fmpz_mod_mat(fmpz_mod_mat_t col, const fq_t a, const fq_ctx_t ctx)
    void fq_set_fmpz_mod_mat(fq_t a, const fmpz_mod_mat_t col, const fq_ctx_t ctx)
    int fq_is_zero(const fq_t op, const fq_ctx_t ctx)
    int fq_is_one(const fq_t op, const fq_ctx_t ctx)
    int fq_equal(const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
    int fq_is_invertible(const fq_t op, const fq_ctx_t ctx)
    int fq_is_invertible_f(fq_t f, const fq_t op, const fq_ctx_t ctx)
    void _fq_trace(fmpz_t rop, const fmpz_struct *op, slong len, const fq_ctx_t ctx)
    void fq_trace(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
    void _fq_norm(fmpz_t rop, const fmpz_struct *op, slong len, const fq_ctx_t ctx)
    void fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
    void _fq_frobenius(fmpz_struct *rop, const fmpz_struct *op, slong len, slong e, const fq_ctx_t ctx)
    void fq_frobenius(fq_t rop, const fq_t op, slong e, const fq_ctx_t ctx)
    int fq_multiplicative_order(fmpz_t ord, const fq_t op, const fq_ctx_t ctx)
    int fq_is_primitive(const fq_t op, const fq_ctx_t ctx)
    void fq_bit_pack(fmpz_t f, const fq_t op, flint_bitcnt_t bit_size, const fq_ctx_t ctx)
    void fq_bit_unpack(fq_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_ctx_t ctx)
