cdef extern from "mpir.h":
    ctypedef unsigned long ulong
    ctypedef unsigned long mp_limb_t
    ctypedef long mp_size_t
    ctypedef long mp_exp_t
    ctypedef unsigned long mp_limb_t
    ctypedef mp_limb_t* mp_ptr
    ctypedef mp_limb_t* mp_srcptr
    ctypedef unsigned long mp_bitcnt_t

ctypedef long fmpz_struct

cdef extern from "flint.h":
    ctypedef void * flint_rand_t
    void flint_randinit(flint_rand_t state)
    void flint_randclear(flint_rand_t state)

cdef extern from "nmod_vec.h":
    ctypedef struct nmod_t:
       mp_limb_t n
       mp_limb_t ninv
       mp_bitcnt_t norm
    void nmod_init(nmod_t * mod, mp_limb_t n)
    mp_limb_t nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_neg(mp_limb_t a, nmod_t mod)
    mp_limb_t nmod_mul(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)

cdef extern from "nmod_poly.h":
    ctypedef struct nmod_poly_struct:
        mp_ptr coeffs
        long alloc
        long length
        nmod_t mod
    ctypedef nmod_poly_struct nmod_poly_t[1]

    ctypedef struct nmod_poly_factor_struct:
        nmod_poly_struct *p
        long *exp
        long num
        long alloc
    ctypedef nmod_poly_factor_struct nmod_poly_factor_t[1]

    void nmod_poly_init(nmod_poly_t poly, mp_limb_t n)
    void nmod_poly_init_preinv(nmod_poly_t poly, mp_limb_t n, mp_limb_t ninv)
    void nmod_poly_init2(nmod_poly_t poly, mp_limb_t n, long alloc)
    void nmod_poly_init2_preinv(nmod_poly_t poly, mp_limb_t n, mp_limb_t ninv, long alloc)
    void nmod_poly_realloc(nmod_poly_t poly, long alloc)
    void nmod_poly_clear(nmod_poly_t poly)
    void nmod_poly_fit_length(nmod_poly_t poly, long alloc)
    long nmod_poly_length(nmod_poly_t poly)
    long nmod_poly_degree(nmod_poly_t poly)
    mp_limb_t nmod_poly_modulus(nmod_poly_t poly)
    mp_bitcnt_t nmod_poly_max_bits(nmod_poly_t poly)
    void nmod_poly_set(nmod_poly_t a, nmod_poly_t b)
    void nmod_poly_swap(nmod_poly_t poly1, nmod_poly_t poly2)
    void nmod_poly_zero(nmod_poly_t res)
    void nmod_poly_truncate(nmod_poly_t poly, long len)
    void nmod_poly_reverse(nmod_poly_t output, nmod_poly_t input, long m)
    void nmod_poly_randtest(nmod_poly_t poly, flint_rand_t state, long len)
    ulong nmod_poly_get_coeff_ui(nmod_poly_t poly, ulong j)
    void nmod_poly_set_coeff_ui(nmod_poly_t poly, ulong j, ulong c)
    char * nmod_poly_get_str(nmod_poly_t poly)
    int nmod_poly_set_str(char * s, nmod_poly_t poly)
    int nmod_poly_print(nmod_poly_t a)
    int nmod_poly_equal(nmod_poly_t a, nmod_poly_t b)
    int nmod_poly_is_zero(nmod_poly_t poly)
    void nmod_poly_shift_left(nmod_poly_t res, nmod_poly_t poly, long k)
    void nmod_poly_shift_right(nmod_poly_t res, nmod_poly_t poly, long k)
    void nmod_poly_add(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    void nmod_poly_sub(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    void nmod_poly_neg(nmod_poly_t res, nmod_poly_t poly1)
    void nmod_poly_scalar_mul_nmod(nmod_poly_t res, nmod_poly_t poly1, mp_limb_t c)
    void nmod_poly_make_monic(nmod_poly_t output, nmod_poly_t input)
    void nmod_poly_mul(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    void nmod_poly_mullow(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2, long trunc)
    void nmod_poly_pow(nmod_poly_t res, nmod_poly_t poly, ulong e)
    void nmod_poly_pow_trunc(nmod_poly_t res, nmod_poly_t poly, ulong e, long trunc)
    void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R, nmod_poly_t A, nmod_poly_t B)
    void nmod_poly_div(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B)
    void nmod_poly_inv_series(nmod_poly_t Qinv, nmod_poly_t Q, long n)
    void nmod_poly_div_series(nmod_poly_t Q, nmod_poly_t A, nmod_poly_t B, long n)
    void nmod_poly_derivative(nmod_poly_t x_prime, nmod_poly_t x)
    void nmod_poly_integral(nmod_poly_t x_int, nmod_poly_t x)
    mp_limb_t nmod_poly_evaluate_nmod(nmod_poly_t poly, mp_limb_t c)
    void nmod_poly_compose(nmod_poly_t res, nmod_poly_t poly1, nmod_poly_t poly2)
    void nmod_poly_gcd(nmod_poly_t G, nmod_poly_t A, nmod_poly_t B)
    void nmod_poly_xgcd(nmod_poly_t G, nmod_poly_t S, nmod_poly_t T, nmod_poly_t A, nmod_poly_t B)
    void nmod_poly_invsqrt_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_sqrt_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_atan_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_tan_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_asin_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_sin_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_cos_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_asinh_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_atanh_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_sinh_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_cosh_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_tanh_series(nmod_poly_t g, nmod_poly_t h, long n)
    void nmod_poly_log_series(nmod_poly_t res, nmod_poly_t f, long n)
    void nmod_poly_exp_series(nmod_poly_t f, nmod_poly_t h, long n)

    int nmod_poly_is_irreducible(nmod_poly_t f)
    mp_limb_t nmod_poly_factor_with_berlekamp(nmod_poly_factor_t result, nmod_poly_t poly)
    mp_limb_t nmod_poly_factor_with_cantor_zassenhaus(nmod_poly_factor_t result, nmod_poly_t poly)
    mp_limb_t nmod_poly_factor(nmod_poly_factor_t result, nmod_poly_t input)
    void nmod_poly_factor_init(nmod_poly_factor_t fac)
    void nmod_poly_factor_clear(nmod_poly_factor_t fac)

cdef extern from "nmod_mat.h":
    ctypedef struct nmod_mat_struct:
        mp_limb_t * entries
        long r
        long c
        mp_limb_t ** rows
        nmod_t mod
    ctypedef nmod_mat_struct nmod_mat_t[1]
    mp_limb_t nmod_mat_entry(nmod_mat_t mat, long i, long j)
    long nmod_mat_nrows(nmod_mat_t mat)
    long nmod_mat_ncols(nmod_mat_t mat)
    void _nmod_mat_set_mod(nmod_mat_t mat, mp_limb_t n)
    void nmod_mat_init(nmod_mat_t mat, long rows, long cols, mp_limb_t n)
    void nmod_mat_init_set(nmod_mat_t mat, nmod_mat_t src)
    void nmod_mat_clear(nmod_mat_t mat)
    void nmod_mat_randtest(nmod_mat_t mat, flint_rand_t state)
    void nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state)
    void nmod_mat_randrank(nmod_mat_t, flint_rand_t state, long rank)
    void nmod_mat_randops(nmod_mat_t mat, long count, flint_rand_t state)
    void nmod_mat_randtril(nmod_mat_t mat, flint_rand_t state, int unit)
    void nmod_mat_randtriu(nmod_mat_t mat, flint_rand_t state, int unit)
    void nmod_mat_print_pretty(nmod_mat_t mat)
    int nmod_mat_equal(nmod_mat_t mat1, nmod_mat_t mat2)
    int nmod_mat_is_zero(nmod_mat_t mat)
    int nmod_mat_is_empty(nmod_mat_t mat)
    int nmod_mat_is_square(nmod_mat_t mat)
    void nmod_mat_zero(nmod_mat_t mat)
    void nmod_mat_set(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_transpose(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_add(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_sub(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_neg(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_scalar_mul(nmod_mat_t B, nmod_mat_t A, mp_limb_t c)
    void nmod_mat_mul(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_mul_classical(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_mul_strassen(nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_addmul(nmod_mat_t D, nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    void nmod_mat_submul(nmod_mat_t D, nmod_mat_t C, nmod_mat_t A, nmod_mat_t B)
    mp_limb_t nmod_mat_det(nmod_mat_t A)
    long nmod_mat_rank(nmod_mat_t A)
    int nmod_mat_inv(nmod_mat_t B, nmod_mat_t A)
    void nmod_mat_solve_tril(nmod_mat_t X, nmod_mat_t L, nmod_mat_t B, int unit)
    void nmod_mat_solve_triu(nmod_mat_t X, nmod_mat_t U, nmod_mat_t B, int unit)
    long nmod_mat_lu(long * P, nmod_mat_t A, int rank_check)
    int nmod_mat_solve(nmod_mat_t X, nmod_mat_t A, nmod_mat_t B)
    int nmod_mat_solve_vec(mp_ptr x, nmod_mat_t A, mp_srcptr b)
    long nmod_mat_rref(long * P, nmod_mat_t A)
    long nmod_mat_nullspace(nmod_mat_t X, nmod_mat_t A)

cdef extern from "fmpz.h":
    ctypedef fmpz_struct fmpz_t[1]
    int COEFF_IS_MPZ(fmpz_struct v)
    void fmpz_init(fmpz_t op)
    void fmpz_clear(fmpz_t op)
    long fmpz_get_si(fmpz_t f)
    ulong fmpz_get_ui(fmpz_t f)
    void fmpz_set_si(fmpz_t f, long val)
    void fmpz_set_ui(fmpz_t f, ulong val)
    #void fmpz_get_mpz(mpz_t x,  fmpz_t f)
    #void fmpz_set_mpz(fmpz_t f,  mpz_t x)
    int fmpz_set_str(fmpz_t f, char * str, int b)
    int fmpz_abs_fits_ui( fmpz_t f)
    void fmpz_zero(fmpz_t f)
    int fmpz_is_zero(fmpz_t f)
    int fmpz_is_one( fmpz_t f)
    int fmpz_is_pm1( fmpz_t f)
    void fmpz_set(fmpz_t f,  fmpz_t g)
    int fmpz_equal(fmpz_t f,  fmpz_t g)
    int fmpz_read(fmpz_t f)
    int fmpz_print(fmpz_t x)
    size_t fmpz_sizeinbase( fmpz_t f, int b)
    char * fmpz_get_str(char * str, int b,  fmpz_t f)
    void fmpz_swap(fmpz_t f, fmpz_t g)
    int fmpz_cmp( fmpz_t f,  fmpz_t g)
    int fmpz_cmp_ui( fmpz_t f, ulong g)
    int fmpz_cmp_si( fmpz_t f, long g)
    int fmpz_cmpabs( fmpz_t f,  fmpz_t g)
    int fmpz_is_even(fmpz_t f)
    int fmpz_is_odd(fmpz_t f)
    mp_size_t fmpz_size(fmpz_t f)
    int fmpz_sgn(fmpz_t f)
    mp_bitcnt_t fmpz_bits(fmpz_t f)
    void fmpz_neg(fmpz_t f1, fmpz_t f2)
    void fmpz_abs(fmpz_t f1,  fmpz_t f2)
    void fmpz_add(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_sub(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_mul_ui(fmpz_t f,  fmpz_t g, ulong x)
    void fmpz_mul_si(fmpz_t f,  fmpz_t g, long x)
    void fmpz_mul(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_mul_2exp(fmpz_t f,  fmpz_t g, ulong exp)
    void fmpz_add_ui(fmpz_t f,  fmpz_t g, ulong x)
    void fmpz_sub_ui(fmpz_t f,  fmpz_t g, ulong x)
    void fmpz_addmul_ui(fmpz_t f,  fmpz_t g, ulong x)
    void fmpz_submul_ui(fmpz_t f,  fmpz_t g, ulong x)
    void fmpz_addmul(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_submul(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_pow_ui(fmpz_t f,  fmpz_t g, ulong exp)
    void fmpz_powm_ui(fmpz_t f,  fmpz_t g, ulong exp,  fmpz_t m)
    void fmpz_powm(fmpz_t f,  fmpz_t g,  fmpz_t e,  fmpz_t m)
    int fmpz_sqrtmod(fmpz_t b,  fmpz_t a,  fmpz_t p)
    void fmpz_sqrt(fmpz_t f,  fmpz_t g)
    void fmpz_sqrtrem(fmpz_t f, fmpz_t r,  fmpz_t g)
    ulong fmpz_fdiv_ui(fmpz_t g, ulong h)
    ulong fmpz_mod_ui(fmpz_t f,  fmpz_t g, ulong h)
    void fmpz_mod(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_gcd(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_lcm(fmpz_t f,  fmpz_t g,  fmpz_t h)
    int fmpz_invmod(fmpz_t f,  fmpz_t g,  fmpz_t h)
    long fmpz_remove(fmpz_t rop,  fmpz_t op,  fmpz_t f)
    void fmpz_divexact(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_divexact_si(fmpz_t f,  fmpz_t g, long h)
    void fmpz_divexact_ui(fmpz_t f,  fmpz_t g, ulong h)
    void fmpz_cdiv_q(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_cdiv_q_si(fmpz_t f,  fmpz_t g, long h)
    void fmpz_cdiv_q_ui(fmpz_t f,  fmpz_t g, ulong h)
    void fmpz_fdiv_qr(fmpz_t f, fmpz_t s,  fmpz_t g,  fmpz_t h)
    void fmpz_fdiv_q(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_fdiv_r(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_fdiv_q_ui(fmpz_t f,  fmpz_t g, ulong h)
    void fmpz_fdiv_q_si(fmpz_t f,  fmpz_t g, long h)
    void fmpz_fdiv_q_2exp(fmpz_t f,  fmpz_t g, ulong exp)
    void fmpz_tdiv_q(fmpz_t f,  fmpz_t g,  fmpz_t h)
    void fmpz_tdiv_q_ui(fmpz_t f,  fmpz_t g, ulong h)
    void fmpz_tdiv_q_si(fmpz_t f,  fmpz_t g, long h)
    double fmpz_get_d_2exp(long * exp,  fmpz_t f)
    void fmpz_mul2_uiui(fmpz_t f,  fmpz_t g, ulong h1, ulong h2)
    void fmpz_divexact2_uiui(fmpz_t f,  fmpz_t g, ulong h1, ulong h2)
    void fmpz_fac_ui(fmpz_t f, ulong n)
    void fmpz_bin_uiui(fmpz_t res, ulong n, ulong k)
    void fmpz_CRT_ui(fmpz_t out,  fmpz_t r1,  fmpz_t m1, ulong r2, ulong m2)
    void fmpz_CRT_ui_unsigned(fmpz_t out,  fmpz_t r1,  fmpz_t m1, ulong r2, ulong m2)
    void fmpz_set_ui_mod(fmpz_t f, mp_limb_t x, mp_limb_t m)

cdef extern from "fmpz_factor.h":
    ctypedef struct fmpz_factor_struct:
        int sign
        fmpz_struct * p
        fmpz_struct * exp
        long alloc
        long num
    ctypedef fmpz_factor_struct fmpz_factor_t[1]
    void fmpz_factor_init(fmpz_factor_t factor)
    void fmpz_factor_clear(fmpz_factor_t factor)
    void fmpz_factor(fmpz_factor_t factor, fmpz_t n)

cdef extern from "fmpz_poly.h":
    ctypedef struct fmpz_poly_struct:
        fmpz_struct * coeffs
        long alloc
        long length
    ctypedef fmpz_poly_struct fmpz_poly_t[1]

    ctypedef struct fmpz_poly_factor_struct:
        fmpz_struct c
        fmpz_poly_struct *p
        long *exp
        long num
        long alloc
    ctypedef fmpz_poly_factor_struct fmpz_poly_factor_t[1]

    void fmpz_poly_init(fmpz_poly_t poly)
    void fmpz_poly_init2(fmpz_poly_t poly, long alloc)
    void fmpz_poly_realloc(fmpz_poly_t poly, long alloc)
    void fmpz_poly_fit_length(fmpz_poly_t poly, long len)
    void fmpz_poly_clear(fmpz_poly_t poly)
    void _fmpz_poly_normalise(fmpz_poly_t poly)
    void _fmpz_poly_set_length(fmpz_poly_t poly, long newlen)
    long fmpz_poly_length(fmpz_poly_t poly)
    long fmpz_poly_degree(fmpz_poly_t poly)
    ulong fmpz_poly_max_limbs(fmpz_poly_t poly)
    long fmpz_poly_max_bits(fmpz_poly_t poly)
    void fmpz_poly_set(fmpz_poly_t poly1, fmpz_poly_t poly2)
    void fmpz_poly_set_ui(fmpz_poly_t poly, ulong c)
    void fmpz_poly_set_si(fmpz_poly_t poly, long c)
    void fmpz_poly_set_fmpz(fmpz_poly_t poly, fmpz_t c)
    #void fmpz_poly_set_mpz(fmpz_poly_t poly, mpz_t c)
    int fmpz_poly_set_str(fmpz_poly_t poly, char * str)
    char * fmpz_poly_get_str(fmpz_poly_t poly)
    char * fmpz_poly_get_str_pretty(fmpz_poly_t poly, char * x)
    void fmpz_poly_zero(fmpz_poly_t poly)
    void fmpz_poly_one(fmpz_poly_t poly)
    void fmpz_poly_zero_coeffs(fmpz_poly_t poly, long i, long j)
    void fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2)
    void fmpz_poly_reverse(fmpz_poly_t res, fmpz_poly_t poly, long n)
    void fmpz_poly_truncate(fmpz_poly_t poly, long newlen)
    #void fmpz_poly_randtest(fmpz_poly_t f, flint_rand_t state, long len, mp_bitcnt_t bits)
    #void fmpz_poly_randtest_unsigned(fmpz_poly_t f, flint_rand_t state, long len, mp_bitcnt_t bits)
    #void fmpz_poly_randtest_not_zero(fmpz_poly_t f, flint_rand_t state, long len, mp_bitcnt_t bits)
    long fmpz_poly_get_coeff_si(fmpz_poly_t poly, long n)
    void fmpz_poly_set_coeff_si(fmpz_poly_t poly, long n, long x)
    ulong fmpz_poly_get_coeff_ui(fmpz_poly_t poly, long n)
    void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, long n, ulong x)
    void fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, long n, fmpz_t x)
    void fmpz_poly_get_coeff_fmpz(fmpz_t x, fmpz_poly_t poly, long n)
    int fmpz_poly_equal(fmpz_poly_t poly1, fmpz_poly_t poly2)
    int fmpz_poly_is_zero(fmpz_poly_t poly)
    int fmpz_poly_is_one(fmpz_poly_t op)
    int fmpz_poly_is_unit(fmpz_poly_t op)
    void fmpz_poly_add(fmpz_poly_t res, fmpz_poly_t poly1, fmpz_poly_t poly2)
    void fmpz_poly_sub(fmpz_poly_t res, fmpz_poly_t poly1, fmpz_poly_t poly2)
    void fmpz_poly_neg(fmpz_poly_t res, fmpz_poly_t poly)
    void fmpz_poly_scalar_mul_ui(fmpz_poly_t poly1, fmpz_poly_t poly2, ulong x)
    void fmpz_poly_scalar_mul_si(fmpz_poly_t poly1, fmpz_poly_t poly2, long x)
    void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t poly1, fmpz_poly_t poly2, fmpz_t x)
    void fmpz_poly_scalar_addmul_fmpz(fmpz_poly_t poly1, fmpz_poly_t poly2, fmpz_t x)
    void fmpz_poly_scalar_submul_fmpz(fmpz_poly_t poly1, fmpz_poly_t poly2, fmpz_t x)
    void fmpz_poly_scalar_fdiv_ui(fmpz_poly_t poly1, fmpz_poly_t poly2, ulong x)
    void fmpz_poly_scalar_fdiv_si(fmpz_poly_t poly1, fmpz_poly_t poly2, long x)
    void fmpz_poly_scalar_fdiv_fmpz(fmpz_poly_t poly1, fmpz_poly_t poly2, fmpz_t x)
    void fmpz_poly_scalar_tdiv_ui(fmpz_poly_t poly1, fmpz_poly_t poly2, ulong x)
    void fmpz_poly_scalar_tdiv_si(fmpz_poly_t poly1, fmpz_poly_t poly2, long x)
    void fmpz_poly_scalar_tdiv_fmpz(fmpz_poly_t poly1, fmpz_poly_t poly2, fmpz_t x)
    void fmpz_poly_scalar_divexact_ui(fmpz_poly_t poly1, fmpz_poly_t poly2, ulong x)
    void fmpz_poly_scalar_divexact_si(fmpz_poly_t poly1, fmpz_poly_t poly2, long x)
    void fmpz_poly_scalar_divexact_fmpz(fmpz_poly_t poly1, fmpz_poly_t poly2, fmpz_t x)
    void fmpz_poly_mul(fmpz_poly_t res, fmpz_poly_t poly1, fmpz_poly_t poly2)
    void fmpz_poly_mullow(fmpz_poly_t res, fmpz_poly_t poly1, fmpz_poly_t poly2, long n)
    void fmpz_poly_mulhigh_n(fmpz_poly_t res, fmpz_poly_t poly1, fmpz_poly_t poly2, long n)
    void fmpz_poly_pow(fmpz_poly_t res, fmpz_poly_t poly, ulong e)
    void fmpz_poly_pow_trunc(fmpz_poly_t res, fmpz_poly_t poly, ulong e, long n)
    void fmpz_poly_shift_left(fmpz_poly_t res, fmpz_poly_t poly, long n)
    void fmpz_poly_shift_right(fmpz_poly_t res, fmpz_poly_t poly, long n)
    void fmpz_poly_2norm(fmpz_t res, fmpz_poly_t poly)
    void fmpz_poly_gcd(fmpz_poly_t res, fmpz_poly_t poly1, fmpz_poly_t poly2)
    void fmpz_poly_xgcd(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t, fmpz_poly_t poly1, fmpz_poly_t poly2)
    void fmpz_poly_content(fmpz_t res, fmpz_poly_t poly)
    void fmpz_poly_primitive_part(fmpz_poly_t res, fmpz_poly_t poly)
    void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B)
    void fmpz_poly_div(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B)
    void fmpz_poly_rem(fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B)
    fmpz_poly_inv_series(fmpz_poly_t Qinv, fmpz_poly_t Q, long n)
    void fmpz_poly_div_series(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B, long n)
    int fmpz_poly_divides(fmpz_poly_t q, fmpz_poly_t a, fmpz_poly_t b)
    void fmpz_poly_derivative(fmpz_poly_t res, fmpz_poly_t poly)
    void fmpz_poly_evaluate_fmpz(fmpz_t res, fmpz_poly_t f, fmpz_t a)
    #void fmpz_poly_evaluate_mpq(mpq_t res, fmpz_poly_t f, mpq_t a)
    mp_limb_t fmpz_poly_evaluate_mod(fmpz_poly_t poly, mp_limb_t a, mp_limb_t n)
    void fmpz_poly_compose(fmpz_poly_t res, fmpz_poly_t poly1, fmpz_poly_t poly2)
    void fmpz_poly_signature(long * r1, long * r2, fmpz_poly_t poly)
    int fmpz_poly_print(fmpz_poly_t poly)
    int fmpz_poly_print_pretty(fmpz_poly_t poly, char * x)
    void fmpz_poly_evaluate_fmpz_vec(fmpz_struct * res, fmpz_poly_t f, fmpz_struct * a, long n)
    void fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly, fmpz_struct * xs, fmpz_struct * ys, long n)
    void fmpz_poly_get_nmod_poly(nmod_poly_t res, fmpz_poly_t poly)
    void fmpz_poly_set_nmod_poly(fmpz_poly_t res, nmod_poly_t poly)
    void fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t res, nmod_poly_t poly)

cdef extern from "fmpz_poly_factor.h":
    void fmpz_poly_factor_init(fmpz_poly_factor_t fac)
    void fmpz_poly_factor_clear(fmpz_poly_factor_t fac)
    void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t fac, fmpz_poly_t G)

cdef extern from "fmpz_mat.h":
    ctypedef struct fmpz_mat_struct:
        fmpz_struct * entries
        long r
        long c
        fmpz_struct ** rows
    ctypedef fmpz_mat_struct fmpz_mat_t[1]
    fmpz_struct * fmpz_mat_entry(fmpz_mat_t mat, long i, long j)
    long fmpz_mat_nrows(fmpz_mat_t mat)
    long fmpz_mat_ncols(fmpz_mat_t mat)
    void fmpz_mat_init(fmpz_mat_t mat, long rows, long cols)
    void fmpz_mat_init_set(fmpz_mat_t mat,  fmpz_mat_t src)
    void fmpz_mat_swap(fmpz_mat_t mat1, fmpz_mat_t mat2)
    void fmpz_mat_set(fmpz_mat_t mat1,  fmpz_mat_t mat2)
    void fmpz_mat_clear(fmpz_mat_t mat)
    int fmpz_mat_equal(fmpz_mat_t mat1, fmpz_mat_t mat2)
    int fmpz_mat_is_zero( fmpz_mat_t mat)
    int fmpz_mat_is_empty( fmpz_mat_t mat)
    int fmpz_mat_is_square( fmpz_mat_t mat)
    void fmpz_mat_zero(fmpz_mat_t mat)
    void fmpz_mat_one(fmpz_mat_t mat)
    void fmpz_mat_get_nmod_mat(nmod_mat_t mat, fmpz_mat_t mat2)
    void fmpz_mat_randbits(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpz_mat_randtest(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpz_mat_randtest_unsigned(fmpz_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpz_mat_randrank(fmpz_mat_t mat, flint_rand_t state, long rank, mp_bitcnt_t bits)
    int fmpz_mat_print_pretty( fmpz_mat_t mat)
    void fmpz_mat_transpose(fmpz_mat_t B,  fmpz_mat_t A)
    void fmpz_mat_add(fmpz_mat_t C,  fmpz_mat_t A,  fmpz_mat_t B)
    void fmpz_mat_sub(fmpz_mat_t C,  fmpz_mat_t A,  fmpz_mat_t B)
    void fmpz_mat_neg(fmpz_mat_t B,  fmpz_mat_t A)
    void fmpz_mat_scalar_mul_fmpz(fmpz_mat_t B,  fmpz_mat_t A,  fmpz_t c)
    void fmpz_mat_scalar_mul_si(fmpz_mat_t B,  fmpz_mat_t A, long c)
    void fmpz_mat_scalar_mul_ui(fmpz_mat_t B,  fmpz_mat_t A, ulong c)
    void fmpz_mat_mul(fmpz_mat_t C,  fmpz_mat_t A,  fmpz_mat_t B)
    void fmpz_mat_det(fmpz_t det,  fmpz_mat_t A)
    long fmpz_mat_rank( fmpz_mat_t A)
    void fmpz_mat_inv(fmpz_mat_t B, fmpz_t den, fmpz_mat_t A)


cdef extern from "fmpq.h":
    ctypedef struct fmpq_struct:
        fmpz_struct num
        fmpz_struct den
    ctypedef fmpq_struct fmpq_t[1]
    fmpz_struct * fmpq_numref(fmpq_t x)
    fmpz_struct * fmpq_denref(fmpq_t x)
    void fmpq_init(fmpq_t x)
    void fmpq_clear(fmpq_t x)
    void fmpq_zero(fmpq_t res)
    void fmpq_one(fmpq_t res)
    int fmpq_equal(fmpq_t x, fmpq_t y)
    int fmpq_sgn(fmpq_t x)
    int fmpq_is_zero(fmpq_t x)
    int fmpq_is_one(fmpq_t x)
    void fmpq_set(fmpq_t dest, fmpq_t src)
    void fmpq_neg(fmpq_t dest, fmpq_t src)
    void fmpq_canonicalise(fmpq_t res)
    int fmpq_is_canonical(fmpq_t x)
    void fmpq_set_si(fmpq_t res, long p, ulong q)
    #void fmpq_set_mpq(fmpq_t dest, mpq_t src)
    #void fmpq_get_mpq(mpq_t dest, fmpq_t src)
    void fmpq_print(fmpq_t x)
    void fmpq_randtest(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_randbits(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_add(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_sub(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_mul(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_mul_fmpz(fmpq_t res, fmpq_t op, fmpz_t x)
    void fmpq_addmul(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_submul(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_inv(fmpq_t dest, fmpq_t src)
    void fmpq_div(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_div_fmpz(fmpq_t res, fmpq_t op, fmpz_t x)
    int fmpq_mod_fmpz(fmpz_t res, fmpq_t x, fmpz_t mod)
    int fmpq_reconstruct_fmpz(fmpq_t res, fmpz_t a, fmpz_t m)
    int fmpq_reconstruct_fmpz_2(fmpq_t res, fmpz_t a, fmpz_t m, fmpz_t N, fmpz_t D)
    mp_bitcnt_t fmpq_height_bits(fmpq_t x)
    void fmpq_height(fmpz_t height, fmpq_t x)
    void fmpq_next_calkin_wilf(fmpq_t res, fmpq_t x)
    void fmpq_next_signed_calkin_wilf(fmpq_t res, fmpq_t x)
    void fmpq_next_minimal(fmpq_t res, fmpq_t x)
    void fmpq_next_signed_minimal(fmpq_t res, fmpq_t x)

cdef extern from "fmpq_poly.h":
    ctypedef struct fmpq_poly_struct:
        fmpz_struct * coeffs
        fmpz_t den
        long alloc
        long length
    ctypedef fmpq_poly_struct fmpq_poly_t[1]
    void fmpq_poly_init(fmpq_poly_t poly)
    void fmpq_poly_init2(fmpq_poly_t poly, long alloc)
    void fmpq_poly_realloc(fmpq_poly_t poly, long alloc)
    void fmpq_poly_fit_length(fmpq_poly_t poly, long len)
    void fmpq_poly_clear(fmpq_poly_t poly)
    void _fmpq_poly_normalise(fmpq_poly_t poly)
    void fmpq_poly_canonicalise(fmpq_poly_t poly)
    int fmpq_poly_is_canonical(fmpq_poly_t poly)
    fmpz_struct * fmpq_poly_numref(fmpq_poly_t poly)
    fmpz_struct * fmpq_poly_denref(fmpq_poly_t poly)
    void fmpq_poly_get_numerator(fmpz_poly_t res, fmpq_poly_t poly)
    long fmpq_poly_degree(fmpq_poly_t poly)
    long fmpq_poly_length(fmpq_poly_t poly)
    void fmpq_poly_randtest(fmpq_poly_t f, flint_rand_t state, long len, mp_bitcnt_t bits_in)
    void fmpq_poly_randtest_unsigned(fmpq_poly_t f, flint_rand_t state, long len, mp_bitcnt_t bits_in)
    void fmpq_poly_randtest_not_zero(fmpq_poly_t f, flint_rand_t state, long len, mp_bitcnt_t bits_in)
    void fmpq_poly_set(fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_set_si(fmpq_poly_t poly, long x)
    void fmpq_poly_set_ui(fmpq_poly_t poly, ulong x)
    void fmpq_poly_set_fmpz(fmpq_poly_t poly, fmpz_t x)
    void fmpq_poly_set_fmpq(fmpq_poly_t poly, fmpq_t x)
    #void fmpq_poly_set_mpz(fmpq_poly_t poly, mpz_t x)
    #void fmpq_poly_set_mpq(fmpq_poly_t poly, mpq_t x)
    void fmpq_poly_set_fmpz_poly(fmpq_poly_t rop, fmpz_poly_t op)
    #void _fmpq_poly_set_array_mpq(fmpz_struct * poly, fmpz_t den, mpq_t * a, long n)
    #void fmpq_poly_set_array_mpq(fmpq_poly_t poly, mpq_t * a, long n)
    int fmpq_poly_set_str(fmpq_poly_t poly, char * str)
    char * fmpq_poly_get_str(fmpq_poly_t poly)
    char * fmpq_poly_get_str_pretty(fmpq_poly_t poly, char * var)
    void fmpq_poly_zero(fmpq_poly_t poly)
    void fmpq_poly_one(fmpq_poly_t poly)
    void fmpq_poly_neg(fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_inv(fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_swap(fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_truncate(fmpq_poly_t poly, long n)
    void fmpq_poly_get_coeff_fmpq(fmpq_t x, fmpq_poly_t poly, long n)
    void fmpq_poly_set_coeff_si(fmpq_poly_t poly, long n, long x)
    void fmpq_poly_set_coeff_ui(fmpq_poly_t poly, long n, ulong x)
    void fmpq_poly_set_coeff_fmpz(fmpq_poly_t poly, long n, fmpz_t x)
    void fmpq_poly_set_coeff_fmpq(fmpq_poly_t poly, long n, fmpq_t x)
    #void fmpq_poly_set_coeff_mpz(fmpq_poly_t poly, long n, mpz_t x)
    #void fmpq_poly_set_coeff_mpq(fmpq_poly_t poly, long n, mpq_t x)
    int fmpq_poly_equal(fmpq_poly_t poly1, fmpq_poly_t poly2)
    int fmpq_poly_cmp(fmpq_poly_t left, fmpq_poly_t right)
    int fmpq_poly_is_zero(fmpq_poly_t poly)
    int fmpq_poly_is_one(fmpq_poly_t poly)
    void fmpq_poly_add(fmpq_poly_t res, fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_sub(fmpq_poly_t res, fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, fmpq_poly_t op, long c)
    void fmpq_poly_scalar_mul_ui(fmpq_poly_t rop, fmpq_poly_t op, ulong c)
    void fmpq_poly_scalar_mul_fmpz(fmpq_poly_t rop, fmpq_poly_t op, fmpz_t c)
    #void fmpq_poly_scalar_mul_mpz(fmpq_poly_t rop, fmpq_poly_t op, mpz_t c)
    #void fmpq_poly_scalar_mul_mpq(fmpq_poly_t rop, fmpq_poly_t op, mpq_t c)
    void fmpq_poly_scalar_div_si(fmpq_poly_t rop, fmpq_poly_t op, long c)
    void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, fmpq_poly_t op, ulong c)
    void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, fmpq_poly_t op, fmpz_t c)
    #void fmpq_poly_scalar_div_mpz(fmpq_poly_t rop, fmpq_poly_t op, mpz_t c)
    #void fmpq_poly_scalar_div_mpq(fmpq_poly_t rop, fmpq_poly_t op, mpq_t c)
    void fmpq_poly_mul(fmpq_poly_t res, fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_mullow(fmpq_poly_t res, fmpq_poly_t poly1, fmpq_poly_t poly2, long n)
    void fmpq_poly_addmul(fmpq_poly_t rop, fmpq_poly_t op1, fmpq_poly_t op2)
    void fmpq_poly_submul(fmpq_poly_t rop, fmpq_poly_t op1, fmpq_poly_t op2)
    void fmpq_poly_pow(fmpq_poly_t rpoly, fmpq_poly_t poly, ulong e)
    void fmpq_poly_shift_left(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_shift_right(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_divrem(fmpq_poly_t Q, fmpq_poly_t R, fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_div(fmpq_poly_t Q, fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_gcd(fmpq_poly_t g, fmpq_poly_t a, fmpq_poly_t b)
    void fmpq_poly_rem(fmpq_poly_t R, fmpq_poly_t poly1, fmpq_poly_t poly2)
    void fmpq_poly_inv_series(fmpq_poly_t Qinv, fmpq_poly_t Q, long n)
    void fmpq_poly_div_series(fmpq_poly_t Q, fmpq_poly_t A, fmpq_poly_t B, long n)
    void fmpq_poly_derivative(fmpq_poly_t res, fmpq_poly_t poly)
    void fmpq_poly_integral(fmpq_poly_t res, fmpq_poly_t poly)
    void fmpq_poly_invsqrt_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_sqrt_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_log_series(fmpq_poly_t res, fmpq_poly_t f, long n)
    void fmpq_poly_exp_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_atan_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_atanh_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_asin_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_asinh_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_tan_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_sin_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_cos_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_sinh_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_cosh_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_tanh_series(fmpq_poly_t res, fmpq_poly_t poly, long n)
    void fmpq_poly_evaluate_fmpz(fmpq_t res, fmpq_poly_t poly, fmpz_t a)
    void fmpq_poly_evaluate_fmpq(fmpq_t res, fmpq_poly_t poly, fmpq_t a)
    # void fmpq_poly_evaluate_mpq(mpq_t res, fmpq_poly_t poly, mpq_t a)
    void fmpq_poly_interpolate_fmpz_vec(fmpq_poly_t poly, fmpz_struct * xs, fmpz_struct * ys, long n)
    void fmpq_poly_compose(fmpq_poly_t res, fmpq_poly_t poly1, fmpq_poly_t poly2)
    #void fmpq_poly_rescale(fmpq_poly_t res, fmpq_poly_t poly, mpq_t x)
    #void fmpq_poly_content(mpq_t res, fmpq_poly_t poly)
    void fmpq_poly_primitive_part(fmpq_poly_t res, fmpq_poly_t poly)
    int fmpq_poly_is_monic(fmpq_poly_t poly)
    void fmpq_poly_make_monic(fmpq_poly_t res, fmpq_poly_t poly)
    int fmpq_poly_is_squarefree(fmpq_poly_t poly)
    int fmpq_poly_debug(fmpq_poly_t poly)
    #int fmpq_poly_fprint(FILE * file, fmpq_poly_t poly)
    int fmpq_poly_print(fmpq_poly_t poly)
    int fmpq_poly_print_pretty(fmpq_poly_t poly, char * var)
    #int fmpq_poly_fread(FILE * file, fmpq_poly_t poly)
    int fmpq_poly_read(fmpq_poly_t poly)

cdef extern from "fmpq_mat.h":
    ctypedef struct fmpq_mat_struct:
        fmpq_struct * entries
        long r
        long c
        fmpq_struct ** rows
    ctypedef fmpq_mat_struct fmpq_mat_t[1]
    fmpq_struct * fmpq_mat_entry(fmpq_mat_t mat, long i, long j)
    fmpz_struct * fmpq_mat_entry_num(fmpq_mat_t mat, long i, long j)
    fmpz_struct * fmpq_mat_entry_den(fmpq_mat_t mat, long i, long j)
    long fmpq_mat_nrows(fmpq_mat_t mat)
    long fmpq_mat_ncols(fmpq_mat_t mat)
    void fmpq_mat_init(fmpq_mat_t mat, long rows, long cols)
    void fmpq_mat_clear(fmpq_mat_t mat)
    void fmpq_mat_print(fmpq_mat_t mat)
    void fmpq_mat_randbits(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_mat_randtest(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_mat_hilbert_matrix(fmpq_mat_t mat)
    void fmpq_mat_set(fmpq_mat_t dest, fmpq_mat_t src)
    void fmpq_mat_zero(fmpq_mat_t mat)
    void fmpq_mat_one(fmpq_mat_t mat)
    void fmpq_mat_add(fmpq_mat_t mat, fmpq_mat_t mat1, fmpq_mat_t mat2)
    void fmpq_mat_sub(fmpq_mat_t mat, fmpq_mat_t mat1, fmpq_mat_t mat2)
    void fmpq_mat_neg(fmpq_mat_t rop, fmpq_mat_t op)
    void fmpq_mat_scalar_mul_fmpz(fmpq_mat_t rop, fmpq_mat_t op, fmpz_t x)
    void fmpq_mat_scalar_div_fmpz(fmpq_mat_t rop, fmpq_mat_t op, fmpz_t x)
    int fmpq_mat_equal(fmpq_mat_t mat1, fmpq_mat_t mat2)
    int fmpq_mat_is_integral(fmpq_mat_t mat)
    int fmpq_mat_is_zero(fmpq_mat_t mat)
    int fmpq_mat_is_empty(fmpq_mat_t mat)
    int fmpq_mat_is_square(fmpq_mat_t mat)
    int fmpq_mat_get_fmpz_mat(fmpz_mat_t dest, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_entrywise(fmpz_mat_t num, fmpz_mat_t den, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t num, fmpz_t den, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_t num, fmpz_struct * den, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_colwise(fmpz_mat_t num, fmpz_struct * den, fmpq_mat_t mat)
    void fmpq_mat_get_fmpz_mat_rowwise_2(fmpz_mat_t num, fmpz_mat_t num2, fmpz_struct * den, fmpq_mat_t mat, fmpq_mat_t mat2)
    void fmpq_mat_get_fmpz_mat_mod_fmpz(fmpz_mat_t dest, fmpq_mat_t mat, fmpz_t mod)
    void fmpq_mat_set_fmpz_mat(fmpq_mat_t dest, fmpz_mat_t src)
    void fmpq_mat_set_fmpz_mat_div_fmpz(fmpq_mat_t X, fmpz_mat_t Xmod, fmpz_t div)
    int fmpq_mat_set_fmpz_mat_mod_fmpz(fmpq_mat_t X, fmpz_mat_t Xmod, fmpz_t mod)
    void fmpq_mat_mul(fmpq_mat_t C, fmpq_mat_t A, fmpq_mat_t B)
    void fmpq_mat_mul_fmpz_mat(fmpq_mat_t C, fmpq_mat_t A, fmpz_mat_t B)
    void fmpq_mat_mul_r_fmpz_mat(fmpq_mat_t C, fmpz_mat_t A, fmpq_mat_t B)
    void fmpq_mat_det(fmpq_t det, fmpq_mat_t mat)
    int fmpq_mat_solve_fraction_free(fmpq_mat_t X, fmpq_mat_t A, fmpq_mat_t B)
    int fmpq_mat_solve_dixon(fmpq_mat_t X, fmpq_mat_t A, fmpq_mat_t B)
    int fmpq_mat_inv(fmpq_mat_t B, fmpq_mat_t A)
    long fmpq_mat_rref(fmpq_mat_t B, fmpq_mat_t A)
    void fmpq_mat_transpose(fmpq_mat_t B, fmpq_mat_t A)


cdef extern from "arith.h":
    void arith_number_of_partitions(fmpz_t res, ulong n)
    int arith_moebius_mu(fmpz_t n)
    void arith_divisor_sigma(fmpz_t v, fmpz_t n, ulong k)
    void arith_euler_phi(fmpz_t v, fmpz_t n)
    void arith_bell_number(fmpz_t v, ulong n)
    void arith_euler_number(fmpz_t v, ulong n)
    void arith_bernoulli_number(fmpq_t v, ulong n)
    void arith_stirling_number_1(fmpz_t s, long n, long k)
    void arith_stirling_number_2(fmpz_t s, long n, long k)
    void arith_harmonic_number(fmpq_t v, ulong n)
    void arith_bernoulli_polynomial(fmpq_poly_t v, ulong n)
    void arith_euler_polynomial(fmpq_poly_t v, ulong n)
    void arith_legendre_polynomial(fmpq_poly_t v, ulong n)
    void arith_chebyshev_t_polynomial(fmpz_poly_t v, ulong n)
    void arith_chebyshev_u_polynomial(fmpz_poly_t v, ulong n)
    void arith_cyclotomic_polynomial(fmpz_poly_t v, ulong n)
