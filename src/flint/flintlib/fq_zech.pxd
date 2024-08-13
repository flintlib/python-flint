from flint.flintlib.flint cimport ulong, slong, flint_bitcnt_t, mp_ptr, flint_rand_t, mp_srcptr, mp_limb_t
from flint.flintlib.fmpz cimport fmpz_t, fmpz_struct
from flint.flintlib.nmod_mat cimport nmod_mat_t
from flint.flintlib.nmod_poly cimport nmod_poly_t, nmod_poly_struct
from flint.flintlib.fq_nmod cimport fq_nmod_t, fq_nmod_ctx_t, fq_nmod_ctx_struct


cdef extern from "flint/fq_zech.h":
    # Type definitions **********************************************/
    ctypedef struct fq_zech_struct:
        mp_limb_t value
    ctypedef fq_zech_struct fq_zech_t[1]

    ctypedef struct fq_zech_ctx_struct:
        mp_limb_t qm1              # q - 1 
        mp_limb_t qm1o2            # (q - 1) / 2 or 1 when p == 2 
        mp_limb_t qm1opm1          # (q - 1) / (p - 1) 
        mp_limb_t p
        double ppre
        mp_limb_t prime_root       # primitive root for prime subfield 
        mp_limb_t * zech_log_table
        mp_limb_t * prime_field_table
        mp_limb_t * eval_table

        fq_nmod_ctx_struct * fq_nmod_ctx
        int owns_fq_nmod_ctx
        int is_conway # whether field was generated using Flint Conway tables (assures primitivity) 
    ctypedef fq_zech_ctx_struct fq_zech_ctx_t[1]

    ctypedef struct fq_zech_mat_struct:
        fq_zech_struct * entries
        slong r
        slong s
        fq_zech_struct ** rows
    ctypedef fq_zech_mat_struct fq_zech_mat_t[1]

    ctypedef struct fq_zech_poly_struct:
        fq_zech_struct * coeffs
        slong alloc
        slong length
    ctypedef fq_zech_poly_struct fq_zech_poly_t[1]

    # Parsed from here **********************************************/
    void fq_zech_ctx_init(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    int _fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    void fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    void fq_zech_ctx_init_random(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var)
    void fq_zech_ctx_init_modulus(fq_zech_ctx_t ctx, nmod_poly_t modulus, const char *var)
    int fq_zech_ctx_init_modulus_check(fq_zech_ctx_t ctx, nmod_poly_t modulus, const char *var)
    void fq_zech_ctx_init_fq_nmod_ctx(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctxn)
    int fq_zech_ctx_init_fq_nmod_ctx_check(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctxn)
    void fq_zech_ctx_clear(fq_zech_ctx_t ctx)
    const nmod_poly_struct* fq_zech_ctx_modulus(const fq_zech_ctx_t ctx)
    long fq_zech_ctx_degree(const fq_zech_ctx_t ctx)
    fmpz_struct * fq_zech_ctx_prime(const fq_zech_ctx_t ctx)
    void fq_zech_ctx_order(fmpz_t f, const fq_zech_ctx_t ctx)
    mp_limb_t fq_zech_ctx_order_ui(const fq_zech_ctx_t ctx)
    # int fq_zech_ctx_fprint(FILE * file, const fq_zech_ctx_t ctx)
    void fq_zech_ctx_print(const fq_zech_ctx_t ctx)
    void fq_zech_ctx_randtest(fq_zech_ctx_t ctx)
    void fq_zech_ctx_randtest_reducible(fq_zech_ctx_t ctx)
    void fq_zech_init(fq_zech_t rop, const fq_zech_ctx_t ctx)
    void fq_zech_init2(fq_zech_t rop, const fq_zech_ctx_t ctx)
    void fq_zech_clear(fq_zech_t rop, const fq_zech_ctx_t ctx)
    void _fq_zech_sparse_reduce(mp_ptr R, slong lenR, const fq_zech_ctx_t ctx)
    void _fq_zech_dense_reduce(mp_ptr R, slong lenR, const fq_zech_ctx_t ctx)
    void _fq_zech_reduce(mp_ptr r, slong lenR, const fq_zech_ctx_t ctx)
    void fq_zech_reduce(fq_zech_t rop, const fq_zech_ctx_t ctx)
    void fq_zech_add(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    void fq_zech_sub(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    void fq_zech_sub_one(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
    void fq_zech_neg(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_mul(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    void fq_zech_mul_fmpz(fq_zech_t rop, const fq_zech_t op, const fmpz_t x, const fq_zech_ctx_t ctx)
    void fq_zech_mul_si(fq_zech_t rop, const fq_zech_t op, slong x, const fq_zech_ctx_t ctx)
    void fq_zech_mul_ui(fq_zech_t rop, const fq_zech_t op, ulong x, const fq_zech_ctx_t ctx)
    void fq_zech_sqr(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_div(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    void _fq_zech_inv(mp_ptr *rop, mp_srcptr *op, slong len, const fq_zech_ctx_t ctx)
    void fq_zech_inv(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_gcdinv(fq_zech_t f, fq_zech_t inv, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void _fq_zech_pow(mp_ptr *rop, mp_srcptr *op, slong len, const fmpz_t e, const fq_zech_ctx_t ctx)
    void fq_zech_pow(fq_zech_t rop, const fq_zech_t op, const fmpz_t e, const fq_zech_ctx_t ctx)
    void fq_zech_pow_ui(fq_zech_t rop, const fq_zech_t op, const ulong e, const fq_zech_ctx_t ctx)
    void fq_zech_sqrt(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
    void fq_zech_pth_root(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
    int fq_zech_is_square(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # int fq_zech_fprint_pretty(FILE *file, const fq_zech_t op, const fq_zech_ctx_t ctx)
    int fq_zech_print_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx)
    # void fq_zech_fprint(FILE * file, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_print(const fq_zech_t op, const fq_zech_ctx_t ctx)
    char * fq_zech_get_str(const fq_zech_t op, const fq_zech_ctx_t ctx)
    char * fq_zech_get_str_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_randtest(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    void fq_zech_randtest_not_zero(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    void fq_zech_randtest_dense(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    void fq_zech_rand(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    void fq_zech_rand_not_zero(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
    void fq_zech_set(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_set_si(fq_zech_t rop, const slong x, const fq_zech_ctx_t ctx)
    void fq_zech_set_ui(fq_zech_t rop, const ulong x, const fq_zech_ctx_t ctx)
    void fq_zech_set_fmpz(fq_zech_t rop, const fmpz_t x, const fq_zech_ctx_t ctx)
    void fq_zech_swap(fq_zech_t op1, fq_zech_t op2, const fq_zech_ctx_t ctx)
    void fq_zech_zero(fq_zech_t rop, const fq_zech_ctx_t ctx)
    void fq_zech_one(fq_zech_t rop, const fq_zech_ctx_t ctx)
    void fq_zech_gen(fq_zech_t rop, const fq_zech_ctx_t ctx)
    int fq_zech_get_fmpz(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_get_fq_nmod(fq_nmod_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_set_fq_nmod(fq_zech_t rop, const fq_nmod_t op, const fq_zech_ctx_t ctx)
    void fq_zech_get_nmod_poly(nmod_poly_t a, const fq_zech_t b, const fq_zech_ctx_t ctx)
    void fq_zech_set_nmod_poly(fq_zech_t a, const nmod_poly_t b, const fq_zech_ctx_t ctx)
    void fq_zech_get_nmod_mat(nmod_mat_t col, const fq_zech_t a, const fq_zech_ctx_t ctx)
    void fq_zech_set_nmod_mat(fq_zech_t a, const nmod_mat_t col, const fq_zech_ctx_t ctx)
    int fq_zech_is_zero(const fq_zech_t op, const fq_zech_ctx_t ctx)
    int fq_zech_is_one(const fq_zech_t op, const fq_zech_ctx_t ctx)
    int fq_zech_equal(const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
    int fq_zech_is_invertible(const fq_zech_t op, const fq_zech_ctx_t ctx)
    int fq_zech_is_invertible_f(fq_zech_t f, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_trace(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_norm(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_frobenius(fq_zech_t rop, const fq_zech_t op, slong e, const fq_zech_ctx_t ctx)
    int fq_zech_multiplicative_order(fmpz_t ord, const fq_zech_t op, const fq_zech_ctx_t ctx)
    int fq_zech_is_primitive(const fq_zech_t op, const fq_zech_ctx_t ctx)
    void fq_zech_bit_pack(fmpz_t f, const fq_zech_t op, flint_bitcnt_t bit_size, const fq_zech_ctx_t ctx)
    void fq_zech_bit_unpack(fq_zech_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_zech_ctx_t ctx)
