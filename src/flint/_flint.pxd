# _flint.pxd
#
# Define the contents of the Python, GMP, Flint and Arb headers.

cdef extern from "Python.h":
    ctypedef void PyObject
    ctypedef void PyTypeObject
    ctypedef long Py_ssize_t
    int PyObject_TypeCheck(object, PyTypeObject*)
    int PyInt_Check(PyObject *o)
    int PyLong_Check(PyObject *o)
    long PyInt_AS_LONG(PyObject *io)
    double PyFloat_AS_DOUBLE(PyObject *io)
    Py_ssize_t PyList_GET_SIZE(PyObject *list)
    long PyLong_AsLongAndOverflow(PyObject *pylong, int *overflow)
    long long PyLong_AsLongLongAndOverflow(PyObject *pylong, int *overflow)

DEF FMPZ_UNKNOWN = 0
DEF FMPZ_REF = 1
DEF FMPZ_TMP = 2

#
# Note: ulong and slong are used throughout Flint/Arb. They are expected to be
# 32 bit unsigned and signed integer types on a 32 bit system and 64 bit on a
# 64 bit system. We denote them as unsigned long and long here which would be
# incorrect on 64 bit Windows but the definition here does not matter because
# their actual sizes will be determined by the values from gmp.h and
# flint/flint.h. Their size in bits (32 or 64) is recorded in the FLINT_BITS
# macro which is defined in flint/flint.h.
#

cdef extern from "gmp.h":
    ctypedef unsigned long ulong
    ctypedef unsigned long mp_limb_t
    ctypedef long mp_size_t
    ctypedef long mp_exp_t
    ctypedef mp_limb_t* mp_ptr
    ctypedef mp_limb_t* mp_srcptr
    ctypedef unsigned long mp_bitcnt_t

cdef extern from "flint/fmpz.h":
    ctypedef long slong
    ctypedef ulong flint_bitcnt_t

ctypedef slong fmpz_struct

cdef extern from "flint/flint.h":
    const int FLINT_BITS
    ctypedef void * flint_rand_t
    void flint_randinit(flint_rand_t state)
    void flint_randclear(flint_rand_t state)
    void flint_set_num_threads(long)
    long flint_get_num_threads()
    void flint_cleanup()

cdef extern from *:
    """
    /* FLINT_BITS is not known until C compile time. We need to check if long
     * or long long matches FLINT_BITS to know which CPython function to call.
     */
    #if FLINT_BITS == 32 && LONG_MAX == 2147483647
    #define pylong_as_slong PyLong_AsLongAndOverflow
    #elif FLINT_BITS == 64 && LLONG_MAX == 9223372036854775807
    #define pylong_as_slong PyLong_AsLongLongAndOverflow
    #else
    #error FLINT_BITS does not match width of long or long long.
    #endif
    """
    slong pylong_as_slong(PyObject *pylong, int *overflow)

cdef extern from "flint/nmod_vec.h":
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

cdef extern from "flint/nmod_poly.h":
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

cdef extern from "flint/nmod_mat.h":
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
    long nmod_mat_rref(nmod_mat_t A)
    long nmod_mat_nullspace(nmod_mat_t X, nmod_mat_t A)

cdef extern from "flint/ulong_extras.h":
    ulong n_gcd(ulong n, ulong k)
    int n_is_prime(ulong n)

cdef extern from "flint/fmpz.h":
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
    void fmpz_one(fmpz_t f)
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
    int fmpz_equal_si(const fmpz_t f, long g)
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
    void fmpz_root(fmpz_t r, const fmpz_t f, long n)
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
    int fmpz_moebius_mu(const fmpz_t f)
    void fmpz_fib_ui(fmpz_t f, ulong n)
    void fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong n)
    void fmpz_rfac_uiui(fmpz_t r, ulong x, ulong n)
    void fmpz_primorial(fmpz_t res, ulong n)
    int fmpz_is_perfect_power(fmpz_t root, const fmpz_t f)
    int fmpz_jacobi(const fmpz_t a, const fmpz_t p)
    int fmpz_is_prime(const fmpz_t n)
    int fmpz_is_probabprime(const fmpz_t n)

cdef extern from "flint/fmpz_factor.h":
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
    int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n, ulong start, ulong num_primes)
    void fmpz_factor_expand(fmpz_t n, const fmpz_factor_t factor)
    void _fmpz_factor_append(fmpz_factor_t factor, const fmpz_t p, ulong exp)
    void fmpz_factor_smooth(fmpz_factor_t factor, fmpz_t n, long bits, int proved)

cdef extern from "flint/fmpz_poly.h":
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
    fmpz_struct * fmpz_poly_get_coeff_ptr(fmpz_poly_t poly, long n)
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

    void fmpz_poly_compose_series(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, long n)
    void fmpz_poly_revert_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, long n)

    void fmpz_poly_signature(long * r1, long * r2, fmpz_poly_t poly)
    int fmpz_poly_print(fmpz_poly_t poly)
    int fmpz_poly_print_pretty(fmpz_poly_t poly, char * x)
    void fmpz_poly_evaluate_fmpz_vec(fmpz_struct * res, fmpz_poly_t f, fmpz_struct * a, long n)
    void fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly, fmpz_struct * xs, fmpz_struct * ys, long n)
    void fmpz_poly_get_nmod_poly(nmod_poly_t res, fmpz_poly_t poly)
    void fmpz_poly_set_nmod_poly(fmpz_poly_t res, nmod_poly_t poly)
    void fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t res, nmod_poly_t poly)
    void fmpz_poly_cyclotomic(fmpz_poly_t, ulong)
    void fmpz_poly_cos_minpoly(fmpz_poly_t, ulong)
    void fmpz_poly_swinnerton_dyer(fmpz_poly_t, ulong)
    int fmpz_poly_sqrt(fmpz_poly_t b, const fmpz_poly_t a)

cdef extern from "flint/fmpz_poly_factor.h":
    void fmpz_poly_factor_init(fmpz_poly_factor_t fac)
    void fmpz_poly_factor_clear(fmpz_poly_factor_t fac)
    void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t fac, fmpz_poly_t G)
    void fmpz_poly_factor(fmpz_poly_factor_t fac, fmpz_poly_t G)
    void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, fmpz_poly_t G)

cdef extern from "flint/fmpz_mat.h":
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
    long fmpz_mat_rank(fmpz_mat_t A)
    long fmpz_mat_rref(fmpz_mat_t R, fmpz_t den, fmpz_mat_t A)
    void fmpz_mat_inv(fmpz_mat_t B, fmpz_t den, fmpz_mat_t A)
    int fmpz_mat_solve(fmpz_mat_t X, fmpz_t den, fmpz_mat_t A, fmpz_mat_t B)
    long fmpz_mat_nullspace(fmpz_mat_t res, fmpz_mat_t mat)
    void fmpz_mat_pow(fmpz_mat_t A, fmpz_mat_t B, ulong e)
    int fmpz_mat_is_hadamard(const fmpz_mat_t A)
    int fmpz_mat_hadamard(fmpz_mat_t A)

    void fmpz_mat_hnf(fmpz_mat_t H, const fmpz_mat_t A)
    void fmpz_mat_hnf_transform(fmpz_mat_t H, fmpz_mat_t U, const  fmpz_mat_t A)
    int fmpz_mat_is_in_hnf(const fmpz_mat_t A)
    void fmpz_mat_snf(fmpz_mat_t S, const fmpz_mat_t A)
    int fmpz_mat_is_in_snf(const fmpz_mat_t A)

    void fmpz_mat_charpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
    void fmpz_mat_minpoly(fmpz_poly_t cp, const fmpz_mat_t mat)

cdef extern from "flint/fmpz_lll.h":
    ctypedef struct fmpz_lll_struct:
        double delta
        double eta
        int rt
        int gt

    ctypedef fmpz_lll_struct fmpz_lll_t[1]

    void fmpz_lll_context_init(fmpz_lll_t fl, double delta, double eta, int rt, int gt)
    void fmpz_lll(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)


cdef extern from "flint/fmpq.h":
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
    void fmpq_harmonic_ui(fmpq_t res, ulong n)
    void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k)

cdef extern from "flint/fmpq_poly.h":
    ctypedef struct fmpq_poly_struct:
        fmpz_struct * coeffs
        long alloc
        long length
        fmpz_t den
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
    void fmpq_poly_scalar_mul_fmpq(fmpq_poly_t rop, fmpq_poly_t op, fmpq_t c)
    void fmpq_poly_scalar_div_si(fmpq_poly_t rop, fmpq_poly_t op, long c)
    void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, fmpq_poly_t op, ulong c)
    void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, fmpq_poly_t op, fmpz_t c)
    #void fmpq_poly_scalar_div_mpz(fmpq_poly_t rop, fmpq_poly_t op, mpz_t c)
    void fmpq_poly_scalar_div_fmpq(fmpq_poly_t rop, fmpq_poly_t op, fmpq_t c)
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
    void fmpq_poly_xgcd(fmpq_poly_t G, fmpq_poly_t S, fmpq_poly_t T, const fmpq_poly_t A, const fmpq_poly_t B)
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
    void fmpq_poly_compose_series(fmpq_poly_t res, fmpq_poly_t poly1, fmpq_poly_t poly2, long n)
    void fmpq_poly_revert_series(fmpq_poly_t res, fmpq_poly_t poly1, long n)

cdef extern from "flint/fmpq_mat.h":
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
    int fmpq_mat_solve_fmpz_mat(fmpq_mat_t X, const fmpz_mat_t A, const fmpz_mat_t B)
    int fmpq_mat_inv(fmpq_mat_t B, fmpq_mat_t A)
    long fmpq_mat_rref(fmpq_mat_t B, fmpq_mat_t A)
    void fmpq_mat_transpose(fmpq_mat_t B, fmpq_mat_t A)

    void fmpq_mat_charpoly(fmpq_poly_t cp, const fmpq_mat_t mat)
    void fmpq_mat_minpoly(fmpq_poly_t cp, const fmpq_mat_t mat)


cdef extern from "flint/arith.h":
    void arith_number_of_partitions(fmpz_t res, ulong n)
    int arith_moebius_mu(fmpz_t n)
    void arith_divisor_sigma(fmpz_t v, ulong k, fmpz_t n)
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

cdef extern from "mag.h":
    ctypedef struct mag_struct:
        fmpz_struct exp
        mp_limb_t man
    ctypedef mag_struct mag_t[1]
    ctypedef mag_struct * mag_ptr
    ctypedef const mag_struct * mag_srcptr

    void mag_init(mag_t x)
    void mag_clear(mag_t x)
    void mag_zero(mag_t x)
    void mag_set(mag_t x, const mag_t y)
    void mag_set_ui_2exp_si(mag_t x, ulong v, long e)
    void mag_hypot(mag_t x, const mag_t y, const mag_t z)

cdef extern from "arf.h":
    ctypedef struct arf_struct:
        fmpz_struct exp
        long size
        mp_limb_t d0
        mp_limb_t d1

    ctypedef arf_struct arf_t[1]
    ctypedef arf_struct * arf_ptr
    ctypedef const arf_struct * arf_srcptr

    ctypedef int arf_rnd_t
    cdef arf_rnd_t ARF_RND_DOWN
    cdef arf_rnd_t ARF_RND_NEAR
    cdef arf_rnd_t ARF_RND_FLOOR
    cdef arf_rnd_t ARF_RND_CEIL
    cdef arf_rnd_t ARF_RND_UP

    void arf_init(arf_t x)
    void arf_clear(arf_t x)
    void arf_zero(arf_t x)
    void arf_pos_inf(arf_t x)
    void arf_neg_inf(arf_t x)
    void arf_nan(arf_t x)
    int arf_is_special(const arf_t x)
    int arf_is_zero(const arf_t x)
    int arf_is_pos_inf(const arf_t x)
    int arf_is_neg_inf(const arf_t x)
    int arf_is_nan(const arf_t x)
    int arf_is_normal(const arf_t x)
    int arf_is_finite(const arf_t x)
    int arf_is_inf(const arf_t x)
    void arf_one(arf_t x)
    int arf_is_one(const arf_t x)
    int arf_sgn(const arf_t x)
    int arf_cmp(const arf_t x, const arf_t y)
    int arf_cmpabs(const arf_t x, const arf_t y)
    void arf_swap(arf_t y, arf_t x)
    void arf_set(arf_t y, const arf_t x)
    void arf_neg(arf_t y, const arf_t x)
    void arf_init_set_ui(arf_t x, ulong v)
    void arf_init_set_si(arf_t x, long v)
    void arf_set_ui(arf_t x, ulong v)
    void arf_set_si(arf_t x, long v)
    int arf_cmpabs_ui(const arf_t x, ulong y)
    void arf_init_set_shallow(arf_t z, const arf_t x)
    void arf_init_neg_shallow(arf_t z, const arf_t x)
    void arf_init_set_mag_shallow(arf_t y, const mag_t x)
    void arf_init_neg_mag_shallow(arf_t z, const mag_t x)
    int arf_cmpabs_mag(const arf_t x, const mag_t y)
    int arf_mag_cmpabs(const mag_t x, const arf_t y)
    void arf_set_fmpz(arf_t y, const fmpz_t x)
    int arf_set_round_ui(arf_t x, ulong v, long prec, arf_rnd_t rnd)
    int arf_set_round_si(arf_t x, long v, long prec, arf_rnd_t rnd)
    int arf_set_round_fmpz(arf_t y, const fmpz_t x, long prec, arf_rnd_t rnd)
    int arf_set_round(arf_t y, const arf_t x, long prec, arf_rnd_t rnd)
    int arf_neg_round(arf_t y, const arf_t x, long prec, arf_rnd_t rnd)
    int arf_equal(const arf_t x, const arf_t y)
    void arf_min(arf_t z, const arf_t a, const arf_t b)
    void arf_max(arf_t z, const arf_t a, const arf_t b)
    void arf_abs(arf_t y, const arf_t x)
    long arf_bits(const arf_t x)
    void arf_bot(fmpz_t e, const arf_t x)
    int arf_is_int(const arf_t x)
    int arf_is_int_2exp_si(const arf_t x, long e)
    int arf_cmp_2exp_si(const arf_t x, long e)
    int arf_cmpabs_2exp_si(const arf_t x, long e)
    void arf_set_si_2exp_si(arf_t x, long man, long exp)
    void arf_set_ui_2exp_si(arf_t x, ulong man, long exp)
    void arf_mul_2exp_si(arf_t y, const arf_t x, long e)
    void arf_mul_2exp_fmpz(arf_t y, const arf_t x, const fmpz_t e)
    int arf_set_round_fmpz_2exp(arf_t y, const fmpz_t x, const fmpz_t exp, long prec, arf_rnd_t rnd)
    void arf_abs_bound_lt_2exp_fmpz(fmpz_t b, const arf_t x)
    void arf_abs_bound_le_2exp_fmpz(fmpz_t b, const arf_t x)
    long arf_abs_bound_lt_2exp_si(const arf_t x)
    void arf_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const arf_t x)
    void arf_get_fmpz(fmpz_t z, const arf_t x, arf_rnd_t rnd)
    long arf_get_si(const arf_t x, arf_rnd_t rnd)
    int arf_get_fmpz_fixed_fmpz(fmpz_t y, const arf_t x, const fmpz_t e)
    int arf_get_fmpz_fixed_si(fmpz_t y, const arf_t x, long e)
    void arf_set_fmpz_2exp(arf_t x, const fmpz_t man, const fmpz_t exp)
    void arf_floor(arf_t z, const arf_t x)
    void arf_ceil(arf_t z, const arf_t x)
    void arf_debug(const arf_t x)
    void arf_print(const arf_t x)
    void arf_printd(const arf_t y, long d)
    void arf_randtest(arf_t x, flint_rand_t state, long bits, long mag_bits)
    void arf_randtest_not_zero(arf_t x, flint_rand_t state, long bits, long mag_bits)
    void arf_randtest_special(arf_t x, flint_rand_t state, long bits, long mag_bits)
    int arf_mul(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)
    int arf_neg_mul(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)
    int arf_mul_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_mul_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_mul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_add(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_add_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_add_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_add_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_add_fmpz_2exp(arf_ptr z, arf_srcptr x, const fmpz_t y, const fmpz_t exp, long prec, arf_rnd_t rnd)
    int arf_sub(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_sub_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_sub_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_sub_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_addmul(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_addmul_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_addmul_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_addmul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_submul(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_submul_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_submul_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_submul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_div(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_div_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_ui_div(arf_ptr z, ulong x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_div_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_si_div(arf_ptr z, long x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_div_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_fmpz_div(arf_ptr z, const fmpz_t x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_fmpz_div_fmpz(arf_ptr z, const fmpz_t x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_sqrt(arf_ptr z, arf_srcptr x, long prec, arf_rnd_t rnd)
    int arf_sqrt_ui(arf_t z, ulong x, long prec, arf_rnd_t rnd)
    int arf_sqrt_fmpz(arf_t z, const fmpz_t x, long prec, arf_rnd_t rnd)
    int arf_rsqrt(arf_ptr z, arf_srcptr x, long prec, arf_rnd_t rnd)
    void arf_get_mag(mag_t y, const arf_t x)
    void arf_get_mag_lower(mag_t y, const arf_t x)
    void arf_set_mag(arf_t y, const mag_t x)
    void mag_init_set_arf(mag_t y, const arf_t x)
    void arf_mag_add_ulp(mag_t z, const mag_t x, const arf_t y, long prec)
    void arf_mag_set_ulp(mag_t z, const arf_t y, long prec)
    void arf_get_fmpq(fmpq_t y, const arf_t x)
    int arf_set_fmpq(arf_t y, const fmpq_t x, long prec, arf_rnd_t rnd)
    double arf_get_d(const arf_t x, arf_rnd_t rnd)
    void arf_set_d(arf_t x, double v)

cdef extern from "arb.h":
    ctypedef struct arb_struct:
        arf_struct mid
        mag_struct rad

    ctypedef arb_struct * arb_ptr
    ctypedef const arb_struct * arb_srcptr
    ctypedef arb_struct arb_t[1]

    arf_ptr arb_midref(const arb_t x)
    mag_ptr arb_radref(const arb_t x)

    void arb_init(arb_t x)
    void arb_clear(arb_t x)

    void arb_init(arb_t x)
    void arb_clear(arb_t x)
    int arb_is_exact(const arb_t x)
    int arb_equal(const arb_t x, const arb_t y)
    int arb_eq(const arb_t x, const arb_t y)
    int arb_ne(const arb_t x, const arb_t y)
    int arb_lt(const arb_t x, const arb_t y)
    int arb_le(const arb_t x, const arb_t y)
    int arb_gt(const arb_t x, const arb_t y)
    int arb_ge(const arb_t x, const arb_t y)
    void arb_zero(arb_t x)
    int arb_is_zero(const arb_t x)
    void arb_pos_inf(arb_t x)
    void arb_neg_inf(arb_t x)
    void arb_zero_pm_inf(arb_t x)
    void arb_indeterminate(arb_t x)
    int arb_is_finite(const arb_t x)
    void arb_set(arb_t x, const arb_t y)
    void arb_swap(arb_t x, arb_t y)
    void arb_set_round(arb_t z, const arb_t x, long prec)
    void arb_trim(arb_t y, const arb_t x)
    void arb_neg(arb_t x, const arb_t y)
    void arb_neg_round(arb_t x, const arb_t y, long prec)
    void arb_abs(arb_t x, const arb_t y)
    void arb_sgn(arb_t x, const arb_t y)
    void arb_set_arf(arb_t x, const arf_t y)
    void arb_set_si(arb_t x, long y)
    void arb_set_ui(arb_t x, ulong y)
    void arb_set_fmpz(arb_t x, const fmpz_t y)
    void arb_set_fmpz_2exp(arb_t x, const fmpz_t y, const fmpz_t exp)
    void arb_set_round_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t exp, long prec)
    void arb_set_round_fmpz(arb_t y, const fmpz_t x, long prec)
    int arb_is_one(const arb_t f)
    void arb_one(arb_t f)
    void arb_print(const arb_t x)
    void arb_printd(const arb_t x, long digits)
    void arb_mul_2exp_si(arb_t y, const arb_t x, long e)
    void arb_mul_2exp_fmpz(arb_t y, const arb_t x, const fmpz_t e)
    int arb_is_int(const arb_t x)
    int arb_contains_zero(const arb_t x)
    int arb_is_nonzero(const arb_t x)
    int arb_is_positive(const arb_t x)
    int arb_is_nonnegative(const arb_t x)
    int arb_is_negative(const arb_t x)
    int arb_is_nonpositive(const arb_t x)
    int arb_contains_negative(const arb_t x)
    int arb_contains_nonpositive(const arb_t x)
    int arb_contains_positive(const arb_t x)
    int arb_contains_nonnegative(const arb_t x)
    int arb_contains_int(const arb_t x)
    void arb_get_mag_lower(mag_t z, const arb_t x)
    void arb_get_mag_lower_nonnegative(mag_t z, const arb_t x)
    void arb_get_mag(mag_t z, const arb_t x)
    void arb_get_abs_ubound_arf(arf_t u, const arb_t x, long prec)
    void arb_get_abs_lbound_arf(arf_t u, const arb_t x, long prec)
    void arb_get_ubound_arf(arf_t u, const arb_t z, long prec)
    void arb_get_lbound_arf(arf_t u, const arb_t z, long prec)
    void arb_nonnegative_part(arb_t u, const arb_t x)
    slong arb_rel_error_bits(const arb_t x)
    slong arb_rel_accuracy_bits(const arb_t x)
    long arb_bits(const arb_t x)
    void arb_randtest_exact(arb_t x, flint_rand_t state, long prec, long mag_bits)
    void arb_randtest_wide(arb_t x, flint_rand_t state, long prec, long mag_bits)
    void arb_randtest_precise(arb_t x, flint_rand_t state, long prec, long mag_bits)
    void arb_randtest(arb_t x, flint_rand_t state, long prec, long mag_bits)
    void arb_randtest_special(arb_t x, flint_rand_t state, long prec, long mag_bits)
    void arb_add_error_arf(arb_t x, const arf_t err)
    void arb_add_error_2exp_si(arb_t x, long err)
    void arb_add_error_2exp_fmpz(arb_t x, const fmpz_t err)
    void arb_add_error(arb_t x, const arb_t error)
    void arb_add_error_mag(arb_t x, const mag_t err)
    int arb_contains_arf(const arb_t x, const arf_t y)
    int arb_contains_fmpq(const arb_t x, const fmpq_t y)
    int arb_contains_fmpz(const arb_t x, const fmpz_t y)
    int arb_contains_si(const arb_t x, long y)
    int arb_overlaps(const arb_t x, const arb_t y)
    int arb_contains(const arb_t x, const arb_t y)
    int arb_contains_interior(const arb_t x, const arb_t y)
    void arb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const arb_t x)
    int arb_get_unique_fmpz(fmpz_t z, const arb_t x)
    void arb_get_fmpz_mid_rad_10exp(fmpz_t mid, fmpz_t rad, fmpz_t exp, const arb_t x, long n)

    int arb_set_str(arb_t res, const char * inp, long prec)

    void arb_floor(arb_t z, const arb_t x, long prec)
    void arb_ceil(arb_t z, const arb_t x, long prec)
    void arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, long prec)
    void arb_get_interval_arf(arf_t a, arf_t b, const arb_t x, long prec)
    void arb_union(arb_t z, const arb_t x, const arb_t y, long prec)
    int arb_intersection(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_min(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_max(arb_t z, const arb_t x, const arb_t y, long prec)

    void arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x, long bits)

    void arb_add(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_add_arf(arb_t z, const arb_t x, const arf_t y, long prec)
    void arb_add_ui(arb_t z, const arb_t x, ulong y, long prec)
    void arb_add_si(arb_t z, const arb_t x, long y, long prec)
    void arb_add_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)
    void arb_add_fmpz_2exp(arb_t z, const arb_t x, const fmpz_t man, const fmpz_t exp, long prec)

    void arb_sub(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_sub_arf(arb_t z, const arb_t x, const arf_t y, long prec)
    void arb_sub_ui(arb_t z, const arb_t x, ulong y, long prec)
    void arb_sub_si(arb_t z, const arb_t x, long y, long prec)
    void arb_sub_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

    void arb_mul(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_mul_arf(arb_t z, const arb_t x, const arf_t y, long prec)
    void arb_mul_si(arb_t z, const arb_t x, long y, long prec)
    void arb_mul_ui(arb_t z, const arb_t x, ulong y, long prec)
    void arb_mul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

    void arb_addmul(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_addmul_arf(arb_t z, const arb_t x, const arf_t y, long prec)
    void arb_addmul_si(arb_t z, const arb_t x, long y, long prec)
    void arb_addmul_ui(arb_t z, const arb_t x, ulong y, long prec)
    void arb_addmul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

    void arb_submul(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_submul_arf(arb_t z, const arb_t x, const arf_t y, long prec)
    void arb_submul_si(arb_t z, const arb_t x, long y, long prec)
    void arb_submul_ui(arb_t z, const arb_t x, ulong y, long prec)
    void arb_submul_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)

    void arb_div(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_div_arf(arb_t z, const arb_t x, const arf_t y, long prec)
    void arb_div_si(arb_t z, const arb_t x, long y, long prec)
    void arb_div_ui(arb_t z, const arb_t x, ulong y, long prec)
    void arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, long prec)
    void arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, long prec)
    void arb_ui_div(arb_t z, ulong x, const arb_t y, long prec)

    void arb_inv(arb_t y, const arb_t x, long prec)
    void arb_set_fmpq(arb_t y, const fmpq_t x, long prec)

    void arb_sqrt(arb_t z, const arb_t x, long prec)
    void arb_sqrt_arf(arb_t z, const arf_t x, long prec)
    void arb_sqrt_fmpz(arb_t z, const fmpz_t x, long prec)
    void arb_sqrt_ui(arb_t z, ulong x, long prec)

    void arb_sqrtpos(arb_t z, const arb_t x, long prec)
    void arb_hypot(arb_t z, const arb_t x, const arb_t y, long prec)

    void arb_rsqrt(arb_t z, const arb_t x, long prec)
    void arb_rsqrt_ui(arb_t z, ulong x, long prec)

    void arb_pow_fmpz_binexp(arb_t y, const arb_t b, const fmpz_t e, long prec)
    void arb_pow_fmpz(arb_t y, const arb_t b, const fmpz_t e, long prec)
    void arb_pow_ui(arb_t y, const arb_t b, ulong e, long prec)
    void arb_ui_pow_ui(arb_t y, ulong b, ulong e, long prec)
    void arb_si_pow_ui(arb_t y, long b, ulong e, long prec)
    void arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, long prec)

    void arb_div_2expm1_ui(arb_t z, const arb_t x, ulong n, long prec)
    void arb_pow(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_root(arb_t z, const arb_t x, ulong k, long prec)
    void arb_log(arb_t z, const arb_t x, long prec)
    void arb_log_arf(arb_t z, const arf_t x, long prec)
    void arb_log_ui(arb_t z, ulong x, long prec)
    void arb_log_fmpz(arb_t z, const fmpz_t x, long prec)
    void arb_log1p(arb_t z, const arb_t x, long prec)
    void arb_exp(arb_t z, const arb_t x, long prec)
    void arb_expm1(arb_t z, const arb_t x, long prec)
    void arb_sin(arb_t s, const arb_t x, long prec)
    void arb_cos(arb_t c, const arb_t x, long prec)
    void arb_sin_cos(arb_t s, arb_t c, const arb_t x, long prec)
    void arb_sin_pi(arb_t s, const arb_t x, long prec)
    void arb_cos_pi(arb_t c, const arb_t x, long prec)
    void arb_sin_cos_pi(arb_t s, arb_t c, const arb_t x, long prec)
    void arb_tan(arb_t y, const arb_t x, long prec)
    void arb_cot(arb_t y, const arb_t x, long prec)
    void arb_sec(arb_t y, const arb_t x, long prec)
    void arb_csc(arb_t y, const arb_t x, long prec)
    void arb_tan_pi(arb_t y, const arb_t x, long prec)
    void arb_cot_pi(arb_t y, const arb_t x, long prec)
    void arb_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, long prec)
    void arb_sin_pi_fmpq(arb_t s, const fmpq_t x, long prec)
    void arb_cos_pi_fmpq(arb_t c, const fmpq_t x, long prec)
    void arb_sinc(arb_t s, const arb_t x, long prec)
    void arb_sinc_pi(arb_t s, const arb_t x, long prec)
    void arb_sinh(arb_t z, const arb_t x, long prec)
    void arb_cosh(arb_t z, const arb_t x, long prec)
    void arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, long prec)
    void arb_tanh(arb_t y, const arb_t x, long prec)
    void arb_coth(arb_t y, const arb_t x, long prec)
    void arb_sech(arb_t y, const arb_t x, long prec)
    void arb_csch(arb_t y, const arb_t x, long prec)
    void arb_atan_arf(arb_t z, const arf_t x, long prec)
    void arb_atan(arb_t z, const arb_t x, long prec)
    void arb_atan2(arb_t z, const arb_t b, const arb_t a, long prec)
    void arb_asin(arb_t z, const arb_t x, long prec)
    void arb_acos(arb_t z, const arb_t x, long prec)
    void arb_atanh(arb_t z, const arb_t x, long prec)
    void arb_asinh(arb_t z, const arb_t x, long prec)
    void arb_acosh(arb_t z, const arb_t x, long prec)
    void arb_fac_ui(arb_t z, ulong n, long prec)
    void arb_bin_ui(arb_t z, const arb_t n, ulong k, long prec)
    void arb_bin_uiui(arb_t z, ulong n, ulong k, long prec)
    void arb_fib_fmpz(arb_t z, const fmpz_t n, long prec)
    void arb_fib_ui(arb_t z, ulong n, long prec)
    void arb_const_pi(arb_t z, long prec)
    void arb_const_sqrt_pi(arb_t z, long prec)
    void arb_const_log_sqrt2pi(arb_t z, long prec)
    void arb_const_log2(arb_t z, long prec)
    void arb_const_log10(arb_t z, long prec)
    void arb_const_euler(arb_t z, long prec)
    void arb_const_catalan(arb_t z, long prec)
    void arb_const_e(arb_t z, long prec)
    void arb_const_khinchin(arb_t z, long prec)
    void arb_const_glaisher(arb_t z, long prec)
    void arb_agm(arb_t z, const arb_t x, const arb_t y, long prec)
    void arb_lgamma(arb_t z, const arb_t x, long prec)
    void arb_rgamma(arb_t z, const arb_t x, long prec)
    void arb_gamma(arb_t z, const arb_t x, long prec)
    void arb_gamma_fmpq(arb_t z, const fmpq_t x, long prec)
    void arb_gamma_fmpz(arb_t z, const fmpz_t x, long prec)
    void arb_digamma(arb_t y, const arb_t x, long prec)
    void arb_zeta(arb_t z, const arb_t s, long prec)
    void arb_zeta_ui(arb_t z, ulong n, long prec)
    void arb_bernoulli_ui(arb_t z, ulong n, long prec)
    void arb_bernoulli_fmpz(arb_t z, const fmpz_t n, long prec)
    void arb_bernoulli_poly_ui(arb_t z, ulong n, const arb_t x, long prec)
    void arb_hurwitz_zeta(arb_t z, const arb_t s, const arb_t a, long prec)

    void arb_bell_fmpz(arb_t z, const fmpz_t n, long prec)

    void arb_partitions_fmpz(arb_t z, const fmpz_t n, long prec)
    void arb_partitions_ui(arb_t z, ulong n, long prec)

    void arb_lambertw(arb_t z, const arb_t x, int flags, long prec)


    void arb_rising_ui_bs(arb_t y, const arb_t x, ulong n, long prec)
    void arb_rising_ui_rs(arb_t y, const arb_t x, ulong n, ulong m, long prec)
    void arb_rising_ui_rec(arb_t y, const arb_t x, ulong n, long prec)
    void arb_rising_ui(arb_t z, const arb_t x, ulong n, long prec)
    void arb_rising_fmpq_ui(arb_t y, const fmpq_t x, ulong n, long prec)
    void arb_rising(arb_t y, const arb_t x, const arb_t n, long prec)

    void arb_rising2_ui_rs(arb_t u, arb_t v, const arb_t x, ulong n, ulong m, long prec)
    void arb_rising2_ui_bs(arb_t u, arb_t v, const arb_t x, ulong n, long prec)
    void arb_rising2_ui(arb_t u, arb_t v, const arb_t x, ulong n, long prec)

    void arb_log_ui_from_prev(arb_t s, ulong k, arb_t log_prev, ulong prev, long prec)

    void arb_const_apery(arb_t s, long prec)

    void arb_zeta_ui_asymp(arb_t x, ulong s, long prec)
    void arb_zeta_ui_borwein_bsplit(arb_t x, ulong s, long prec)
    void arb_zeta_ui_euler_product(arb_t z, ulong s, long prec)
    void arb_zeta_ui_bernoulli(arb_t x, ulong n, long prec)
    void arb_zeta_ui_vec_borwein(arb_ptr z, ulong start, long num, ulong step, long prec)
    void arb_zeta_ui(arb_t x, ulong n, long prec)
    void arb_zeta_ui_vec_even(arb_ptr x, ulong start, long num, long prec)
    void arb_zeta_ui_vec_odd(arb_ptr x, ulong start, long num, long prec)
    void arb_zeta_ui_vec(arb_ptr x, ulong start, long num, long prec)
    void arb_bernoulli_ui(arb_t b, ulong n, long prec)
    void arb_bernoulli_ui_zeta(arb_t b, ulong n, long prec)

    void arb_polylog(arb_t w, const arb_t s, const arb_t z, long prec)
    void arb_polylog_si(arb_t w, long s, const arb_t z, long prec)

    void arb_chebyshev_t_ui(arb_t a, ulong n, const arb_t x, long prec)
    void arb_chebyshev_t2_ui(arb_t a, arb_t b, ulong n, const arb_t x, long prec)
    void arb_chebyshev_u_ui(arb_t a, ulong n, const arb_t x, long prec)
    void arb_chebyshev_u2_ui(arb_t a, arb_t b, ulong n, const arb_t x, long prec)

    void arb_root_ui(arb_t z, const arb_t x, ulong k, long prec)

    cdef ulong ARB_STR_MORE
    cdef ulong ARB_STR_NO_RADIUS
    cdef ulong ARB_STR_CONDENSE
    char * arb_get_str(const arb_t x, long n, ulong flags)

cdef extern from "acb.h":
    ctypedef struct acb_struct:
        arb_struct real
        arb_struct imag

    ctypedef acb_struct * acb_ptr
    ctypedef const acb_struct * acb_srcptr
    ctypedef acb_struct acb_t[1]

    arb_ptr acb_realref(const acb_t x)
    arb_ptr acb_imagref(const acb_t x)

    acb_ptr _acb_vec_init(long n)
    void _acb_vec_clear(acb_ptr v, long n)
    void _acb_vec_sort_pretty(acb_ptr vec, long len)
    void acb_printd(const acb_t z, long digits)

    void acb_init(acb_t x)
    void acb_clear(acb_t x)
    int acb_is_zero(const acb_t z)
    int acb_is_one(const acb_t z)
    int acb_is_exact(const acb_t z)
    int acb_is_finite(const acb_t x)
    void acb_indeterminate(acb_t x)
    void acb_zero(acb_t z)
    void acb_one(acb_t z)
    void acb_onei(acb_t z)
    void acb_set(acb_t z, const acb_t x)
    void acb_set_round(acb_t z, const acb_t x, long prec)
    void acb_neg_round(acb_t z, const acb_t x, long prec)
    void acb_swap(acb_t z, acb_t x)
    int acb_equal(const acb_t x, const acb_t y)
    int acb_eq(const acb_t x, const acb_t y)
    int acb_ne(const acb_t x, const acb_t y)
    int acb_overlaps(const acb_t x, const acb_t y)
    int acb_contains_zero(const acb_t x)
    int acb_contains_fmpq(const acb_t x, const fmpq_t y)
    int acb_contains_fmpz(const acb_t x, const fmpz_t y)
    int acb_contains(const acb_t x, const acb_t y)
    int acb_contains_interior(const acb_t x, const acb_t y)
    int acb_get_unique_fmpz(fmpz_t z, const acb_t x)
    int acb_contains_int(const acb_t x)
    void acb_union(acb_t z, const acb_t x, const acb_t y, long prec)
    void acb_set_ui(acb_t z, ulong c)
    void acb_set_si(acb_t z, long c)
    void acb_set_fmpz(acb_t z, const fmpz_t c)
    void acb_set_round_fmpz(acb_t z, const fmpz_t y, long prec)
    void acb_set_fmpq(acb_t z, const fmpq_t c, long prec)
    void acb_set_arb(acb_t z, const arb_t c)
    void acb_set_round_arb(acb_t z, const arb_t x, long prec)
    void acb_trim(acb_t z, const acb_t x)
    void acb_add_error_arf(acb_t x, const arf_t err)
    void acb_add_error_mag(acb_t x, const mag_t err)
    void acb_get_mag(mag_t z, const acb_t x)
    void acb_get_mag_lower(mag_t z, const acb_t x)
    void acb_get_abs_ubound_arf(arf_t u, const acb_t z, long prec)
    void acb_get_abs_lbound_arf(arf_t u, const acb_t z, long prec)
    void acb_get_rad_ubound_arf(arf_t u, const acb_t z, long prec)
    void acb_arg(arb_t r, const acb_t z, long prec)
    void acb_add(acb_t z, const acb_t x, const acb_t y, long prec)
    void acb_sub(acb_t z, const acb_t x, const acb_t y, long prec)
    void acb_add_ui(acb_t z, const acb_t x, ulong c, long prec)
    void acb_sub_ui(acb_t z, const acb_t x, ulong c, long prec)
    void acb_add_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
    void acb_add_arb(acb_t z, const acb_t x, const arb_t y, long prec)
    void acb_sub_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
    void acb_sub_arb(acb_t z, const acb_t x, const arb_t y, long prec)
    void acb_neg(acb_t z, const acb_t x)
    void acb_conj(acb_t z, const acb_t x)
    void acb_abs(arb_t u, const acb_t z, long prec)
    void acb_sgn(acb_t u, const acb_t z, long prec)
    void acb_csgn(arb_t u, const acb_t z)
    void acb_get_real(arb_t u, const acb_t z)
    void acb_get_imag(arb_t u, const acb_t z)

    void acb_real_abs(acb_t res, const acb_t z, int analytic, long prec)
    void acb_real_sgn(acb_t res, const acb_t z, int analytic, long prec)
    void acb_real_heaviside(acb_t res, const acb_t z, int analytic, long prec)
    void acb_real_floor(acb_t res, const acb_t z, int analytic, long prec)
    void acb_real_ceil(acb_t res, const acb_t z, int analytic, long prec)
    void acb_real_max(acb_t res, const acb_t x, const acb_t y, int analytic, long prec)
    void acb_real_min(acb_t res, const acb_t x, const acb_t y, int analytic, long prec)
    void acb_real_sqrtpos(acb_t res, const acb_t z, int analytic, long prec)

    void acb_sqrt_analytic(acb_t res, const acb_t z, int analytic, long prec)
    void acb_rsqrt_analytic(acb_t res, const acb_t z, int analytic, long prec)
    void acb_log_analytic(acb_t res, const acb_t z, int analytic, long prec)
    void acb_pow_analytic(acb_t res, const acb_t z, const acb_t w, int analytic, long prec)

    void acb_mul_ui(acb_t z, const acb_t x, ulong y, long prec)
    void acb_mul_si(acb_t z, const acb_t x, long y, long prec)
    void acb_mul_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
    void acb_mul_arb(acb_t z, const acb_t x, const arb_t y, long prec)
    void acb_mul_onei(acb_t z, const acb_t x)
    void acb_mul(acb_t z, const acb_t x, const acb_t y, long prec)
    void acb_mul_2exp_si(acb_t z, const acb_t x, long e)
    void acb_mul_2exp_fmpz(acb_t z, const acb_t x, const fmpz_t c)
    void acb_addmul(acb_t z, const acb_t x, const acb_t y, long prec)
    void acb_submul(acb_t z, const acb_t x, const acb_t y, long prec)
    void acb_addmul_ui(acb_t z, const acb_t x, ulong y, long prec)
    void acb_addmul_si(acb_t z, const acb_t x, long y, long prec)
    void acb_submul_ui(acb_t z, const acb_t x, ulong y, long prec)
    void acb_submul_si(acb_t z, const acb_t x, long y, long prec)
    void acb_addmul_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
    void acb_submul_fmpz(acb_t z, const acb_t x, const fmpz_t y, long prec)
    void acb_addmul_arb(acb_t z, const acb_t x, const arb_t y, long prec)
    void acb_submul_arb(acb_t z, const acb_t x, const arb_t y, long prec)
    void acb_inv(acb_t z, const acb_t x, long prec)
    void acb_div(acb_t z, const acb_t x, const acb_t y, long prec)
    void acb_div_ui(acb_t z, const acb_t x, ulong c, long prec)
    void acb_div_si(acb_t z, const acb_t x, long c, long prec)
    void acb_div_arb(acb_t z, const acb_t x, const arb_t c, long prec)
    void acb_div_fmpz(acb_t z, const acb_t x, const fmpz_t c, long prec)
    void acb_cube(acb_t y, const acb_t x, long prec)
    void acb_pow_fmpz(acb_t y, const acb_t b, const fmpz_t e, long prec)
    void acb_pow_ui(acb_t y, const acb_t b, ulong e, long prec)
    void acb_pow_si(acb_t y, const acb_t b, long e, long prec)
    void acb_const_pi(acb_t x, long prec)
    void acb_log(acb_t r, const acb_t z, long prec)
    void acb_exp(acb_t r, const acb_t z, long prec)
    void acb_exp_pi_i(acb_t r, const acb_t z, long prec)
    void acb_sin(acb_t r, const acb_t z, long prec)
    void acb_cos(acb_t r, const acb_t z, long prec)
    void acb_sin_cos(acb_t s, acb_t c, const acb_t z, long prec)
    void acb_tan(acb_t r, const acb_t z, long prec)
    void acb_cot(acb_t r, const acb_t z, long prec)
    void acb_sec(acb_t r, const acb_t z, long prec)
    void acb_csc(acb_t r, const acb_t z, long prec)
    void acb_sin_pi(acb_t r, const acb_t z, long prec)
    void acb_cos_pi(acb_t r, const acb_t z, long prec)
    void acb_sin_cos_pi(acb_t s, acb_t c, const acb_t z, long prec)
    void acb_tan_pi(acb_t r, const acb_t z, long prec)
    void acb_cot_pi(acb_t r, const acb_t z, long prec)
    void acb_sinh(acb_t r, const acb_t z, long prec)
    void acb_cosh(acb_t r, const acb_t z, long prec)
    void acb_sinh_cosh(acb_t s, acb_t c, const acb_t z, long prec)
    void acb_tanh(acb_t r, const acb_t z, long prec)
    void acb_coth(acb_t r, const acb_t z, long prec)
    void acb_sech(acb_t r, const acb_t z, long prec)
    void acb_csch(acb_t r, const acb_t z, long prec)
    void acb_sinc(acb_t r, const acb_t z, long prec)
    void acb_sinc_pi(acb_t r, const acb_t z, long prec)
    void acb_pow_arb(acb_t z, const acb_t x, const arb_t y, long prec)
    void acb_pow(acb_t r, const acb_t x, const acb_t y, long prec)
    void acb_sqrt(acb_t y, const acb_t x, long prec)
    void acb_rsqrt(acb_t y, const acb_t x, long prec)
    void acb_rising_ui_bs(acb_t y, const acb_t x, ulong n, long prec)
    void acb_rising_ui_rs(acb_t y, const acb_t x, ulong n, ulong m, long prec)
    void acb_rising_ui_rec(acb_t y, const acb_t x, ulong n, long prec)
    void acb_rising_ui(acb_t z, const acb_t x, ulong n, long prec)
    void acb_rising2_ui_bs(acb_t u, acb_t v, const acb_t x, ulong n, long prec)
    void acb_rising2_ui_rs(acb_t u, acb_t v, const acb_t x, ulong n, ulong m, long prec)
    void acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, long prec)
    void acb_rising_ui_get_mag(mag_t bound, const acb_t s, ulong n)
    void acb_rising(acb_t y, const acb_t x, const acb_t n, long prec)

    void acb_gamma(acb_t y, const acb_t x, long prec)
    void acb_rgamma(acb_t y, const acb_t x, long prec)
    void acb_lgamma(acb_t y, const acb_t x, long prec)
    void acb_digamma(acb_t y, const acb_t x, long prec)
    void acb_zeta(acb_t z, const acb_t s, long prec)
    void acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, long prec)
    void acb_polylog(acb_t w, const acb_t s, const acb_t z, long prec)
    void acb_polylog_si(acb_t w, long s, const acb_t z, long prec)
    void acb_agm1(acb_t m, const acb_t z, long prec)
    void acb_agm1_cpx(acb_ptr m, const acb_t z, long len, long prec)
    void acb_agm(acb_t res, const acb_t a, const acb_t b, long prec)
    void acb_expm1(acb_t r, const acb_t z, long prec)
    void acb_log1p(acb_t r, const acb_t z, long prec)
    void acb_asin(acb_t r, const acb_t z, long prec)
    void acb_acos(acb_t r, const acb_t z, long prec)
    void acb_atan(acb_t r, const acb_t z, long prec)
    void acb_asinh(acb_t r, const acb_t z, long prec)
    void acb_acosh(acb_t r, const acb_t z, long prec)
    void acb_atanh(acb_t r, const acb_t z, long prec)
    void acb_log_sin_pi(acb_t res, const acb_t z, long prec)

    void acb_polygamma(acb_t w, const acb_t s, const acb_t z, long prec)
    void acb_log_barnes_g(acb_t w, const acb_t z, long prec)
    void acb_barnes_g(acb_t w, const acb_t z, long prec)

    void acb_bernoulli_poly_ui(acb_t res, ulong n, const acb_t x, long prec)

    void acb_lambertw(acb_t z, const acb_t x, const fmpz_t k, int flags, long prec)

    long acb_rel_error_bits(const acb_t x)
    long acb_rel_accuracy_bits(const acb_t x)
    long acb_bits(const acb_t x)

    void acb_root_ui(acb_t z, const acb_t x, ulong k, long prec)

cdef extern from "partitions.h":
    void partitions_fmpz_fmpz(fmpz_t, const fmpz_t, int)

cdef extern from "bernoulli.h":
    void bernoulli_fmpq_ui(fmpq_t, ulong)
    void bernoulli_cache_compute(long n)


cdef extern from "arb_poly.h":
    ctypedef struct arb_poly_struct:
        arb_ptr coeffs
        long length
        long alloc

    ctypedef arb_poly_struct arb_poly_t[1]

    void arb_poly_init(arb_poly_t poly)
    void arb_poly_init2(arb_poly_t poly, long len)
    void arb_poly_clear(arb_poly_t poly)
    void arb_poly_fit_length(arb_poly_t poly, long len)
    void _arb_poly_set_length(arb_poly_t poly, long len)
    void _arb_poly_normalise(arb_poly_t poly)
    void arb_poly_swap(arb_poly_t poly1, arb_poly_t poly2)
    void arb_poly_set(arb_poly_t poly, const arb_poly_t src)
    void arb_poly_set_round(arb_poly_t dest, const arb_poly_t src, long prec)
    long arb_poly_length(const arb_poly_t poly)
    long arb_poly_degree(const arb_poly_t poly)
    void arb_poly_zero(arb_poly_t poly)
    void arb_poly_one(arb_poly_t poly)
    void arb_poly_set_coeff_si(arb_poly_t poly, long n, long x)
    void arb_poly_set_coeff_arb(arb_poly_t poly, long n, const arb_t x)
    void arb_poly_get_coeff_arb(arb_t x, const arb_poly_t poly, long n)
    arb_ptr arb_poly_get_coeff_ptr(arb_poly_t poly, long n)
    void _arb_poly_reverse(arb_ptr res, arb_srcptr poly, long len, long n)
    void _arb_poly_shift_right(arb_ptr res, arb_srcptr poly, long len, long n)
    void arb_poly_shift_right(arb_poly_t res, const arb_poly_t poly, long n)
    void _arb_poly_shift_left(arb_ptr res, arb_srcptr poly, long len, long n)
    void arb_poly_shift_left(arb_poly_t res, const arb_poly_t poly, long n)
    void arb_poly_truncate(arb_poly_t poly, long newlen)
    void arb_poly_set_fmpz_poly(arb_poly_t poly, const fmpz_poly_t src, long prec)
    void arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src, long prec)
    void arb_poly_set_arb(arb_poly_t poly, const arb_t c)
    void arb_poly_set_si(arb_poly_t poly, long c)
    int arb_poly_contains(const arb_poly_t poly1, const arb_poly_t poly2)
    int arb_poly_contains_fmpz_poly(const arb_poly_t poly1, const fmpz_poly_t poly2)
    int arb_poly_contains_fmpq_poly(const arb_poly_t poly1, const fmpq_poly_t poly2)
    int arb_poly_equal(const arb_poly_t A, const arb_poly_t B)
    int _arb_poly_overlaps(arb_srcptr poly1, long len1, arb_srcptr poly2, long len2)
    int arb_poly_overlaps(const arb_poly_t poly1, const arb_poly_t poly2)
    void arb_poly_printd(const arb_poly_t poly, long digits)
    void arb_poly_randtest(arb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits)
    void _arb_poly_add(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long prec)
    void arb_poly_add(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)
    void _arb_poly_sub(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long prec)
    void arb_poly_sub(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)
    void arb_poly_neg(arb_poly_t res, const arb_poly_t poly)
    void arb_poly_scalar_mul_2exp_si(arb_poly_t res, const arb_poly_t poly, long c)
    void _arb_poly_mullow(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long n, long prec)
    void arb_poly_mullow(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long len, long prec)
    void _arb_poly_mul(arb_ptr C, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)
    void arb_poly_mul(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)
    void _arb_poly_inv_series(arb_ptr Qinv, arb_srcptr Q, long Qlen, long len, long prec)
    void arb_poly_inv_series(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec)
    void  _arb_poly_div_series(arb_ptr Q, arb_srcptr A, long Alen, arb_srcptr B, long Blen, long n, long prec)
    void arb_poly_div_series(arb_poly_t Q, const arb_poly_t A, const arb_poly_t B, long n, long prec)
    void _arb_poly_div(arb_ptr Q, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)
    void _arb_poly_divrem(arb_ptr Q, arb_ptr R, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)
    void _arb_poly_rem(arb_ptr R, arb_srcptr A, long lenA, arb_srcptr B, long lenB, long prec)
    int arb_poly_divrem(arb_poly_t Q, arb_poly_t R, const arb_poly_t A, const arb_poly_t B, long prec)
    void _arb_poly_div_root(arb_ptr Q, arb_t R, arb_srcptr A, long len, const arb_t c, long prec)
    void _arb_poly_product_roots(arb_ptr poly, arb_srcptr xs, long n, long prec)
    void arb_poly_product_roots(arb_poly_t poly, arb_srcptr xs, long n, long prec)
    arb_ptr * _arb_poly_tree_alloc(long len)
    void _arb_poly_tree_free(arb_ptr * tree, long len)
    void _arb_poly_tree_build(arb_ptr * tree, arb_srcptr roots, long len, long prec)
    void _arb_poly_compose(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long prec)
    void arb_poly_compose(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)
    void _arb_poly_compose_horner(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long prec)
    void arb_poly_compose_horner(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)
    void _arb_poly_compose_divconquer(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long prec)
    void arb_poly_compose_divconquer(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long prec)
    void _arb_poly_compose_series_horner(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long n, long prec)
    void arb_poly_compose_series_horner(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long n, long prec)
    void _arb_poly_compose_series(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long n, long prec)
    void arb_poly_compose_series(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long n, long prec)
    void _arb_poly_revert_series_lagrange(arb_ptr Qinv, arb_srcptr Q, long Qlen, long n, long prec)
    void arb_poly_revert_series_lagrange(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec)
    void _arb_poly_revert_series_newton(arb_ptr Qinv, arb_srcptr Q, long Qlen, long n, long prec)
    void arb_poly_revert_series_newton(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec)
    void _arb_poly_revert_series_lagrange_fast(arb_ptr Qinv, arb_srcptr Q, long Qlen, long n, long prec)
    void arb_poly_revert_series_lagrange_fast(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec)
    void _arb_poly_revert_series(arb_ptr Qinv, arb_srcptr Q, long Qlen, long n, long prec)
    void arb_poly_revert_series(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec)

    void _arb_poly_evaluate_horner(arb_t res, arb_srcptr f, long len, const arb_t a, long prec)
    void arb_poly_evaluate_horner(arb_t res, const arb_poly_t f, const arb_t a, long prec)
    void _arb_poly_evaluate_rectangular(arb_t y, arb_srcptr poly, long len, const arb_t x, long prec)
    void arb_poly_evaluate_rectangular(arb_t res, const arb_poly_t f, const arb_t a, long prec)
    void _arb_poly_evaluate(arb_t res, arb_srcptr f, long len, const arb_t a, long prec)
    void arb_poly_evaluate(arb_t res, const arb_poly_t f, const arb_t a, long prec)
    void _arb_poly_evaluate2_horner(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec)
    void arb_poly_evaluate2_horner(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec)
    void _arb_poly_evaluate2_rectangular(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec)
    void arb_poly_evaluate2_rectangular(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec)
    void _arb_poly_evaluate2(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec)
    void arb_poly_evaluate2(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec)

    void _arb_poly_evaluate_vec_iter(arb_ptr ys, arb_srcptr poly, long plen, arb_srcptr xs, long n, long prec)
    void arb_poly_evaluate_vec_iter(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, long n, long prec)
    void _arb_poly_evaluate_vec_fast_precomp(arb_ptr vs, arb_srcptr poly, long plen, arb_ptr * tree, long len, long prec)
    void _arb_poly_evaluate_vec_fast(arb_ptr ys, arb_srcptr poly, long plen, arb_srcptr xs, long n, long prec)
    void arb_poly_evaluate_vec_fast(arb_ptr ys, const arb_poly_t poly, arb_srcptr xs, long n, long prec)
    void _arb_poly_interpolate_newton(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)
    void arb_poly_interpolate_newton(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)
    void _arb_poly_interpolate_barycentric(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)
    void arb_poly_interpolate_barycentric(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)
    void _arb_poly_interpolation_weights(arb_ptr w, arb_ptr * tree, long len, long prec)
    void _arb_poly_interpolate_fast_precomp(arb_ptr poly, arb_srcptr ys, arb_ptr * tree, arb_srcptr weights, long len, long prec)
    void _arb_poly_interpolate_fast(arb_ptr poly, arb_srcptr xs, arb_srcptr ys, long len, long prec)
    void arb_poly_interpolate_fast(arb_poly_t poly, arb_srcptr xs, arb_srcptr ys, long n, long prec)

    void _arb_poly_derivative(arb_ptr res, arb_srcptr poly, long len, long prec)
    void arb_poly_derivative(arb_poly_t res, const arb_poly_t poly, long prec)
    void _arb_poly_integral(arb_ptr res, arb_srcptr poly, long len, long prec)
    void arb_poly_integral(arb_poly_t res, const arb_poly_t poly, long prec)

    void arb_poly_borel_transform(arb_poly_t res, const arb_poly_t poly, long prec)
    void _arb_poly_borel_transform(arb_ptr res, arb_srcptr poly, long len, long prec)
    void arb_poly_inv_borel_transform(arb_poly_t res, const arb_poly_t poly, long prec)
    void _arb_poly_inv_borel_transform(arb_ptr res, arb_srcptr poly, long len, long prec)
    void _arb_poly_binomial_transform_basecase(arb_ptr b, arb_srcptr a, long alen, long len, long prec)
    void arb_poly_binomial_transform_basecase(arb_poly_t b, const arb_poly_t a, long len, long prec)
    void _arb_poly_binomial_transform_convolution(arb_ptr b, arb_srcptr a, long alen, long len, long prec)
    void arb_poly_binomial_transform_convolution(arb_poly_t b, const arb_poly_t a, long len, long prec)
    void _arb_poly_binomial_transform(arb_ptr b, arb_srcptr a, long alen, long len, long prec)
    void arb_poly_binomial_transform(arb_poly_t b, const arb_poly_t a, long len, long prec)

    void _arb_poly_pow_ui_trunc_binexp(arb_ptr res, arb_srcptr f, long flen, ulong exp, long len, long prec)
    void arb_poly_pow_ui_trunc_binexp(arb_poly_t res, const arb_poly_t poly, ulong exp, long len, long prec)
    void _arb_poly_pow_ui(arb_ptr res, arb_srcptr f, long flen, ulong exp, long prec)
    void arb_poly_pow_ui(arb_poly_t res, const arb_poly_t poly, ulong exp, long prec)
    void _arb_poly_pow_series(arb_ptr h, arb_srcptr f, long flen, arb_srcptr g, long glen, long len, long prec)
    void arb_poly_pow_series(arb_poly_t h, const arb_poly_t f, const arb_poly_t g, long len, long prec)
    void _arb_poly_pow_arb_series(arb_ptr h, arb_srcptr f, long flen, const arb_t g, long len, long prec)
    void arb_poly_pow_arb_series(arb_poly_t h, const arb_poly_t f, const arb_t g, long len, long prec)
    void _arb_poly_rsqrt_series(arb_ptr g, arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_rsqrt_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    void _arb_poly_sqrt_series(arb_ptr g, arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_sqrt_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    void _arb_poly_log_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)
    void arb_poly_log_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    void _arb_poly_atan_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)
    void arb_poly_atan_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    void _arb_poly_asin_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)
    void arb_poly_asin_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    void _arb_poly_acos_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec)
    void arb_poly_acos_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    void _arb_poly_exp_series_basecase(arb_ptr f, arb_srcptr h, long hlen, long n, long prec)
    void arb_poly_exp_series_basecase(arb_poly_t f, const arb_poly_t h, long n, long prec)
    void _arb_poly_exp_series(arb_ptr f, arb_srcptr h, long hlen, long n, long prec)
    void arb_poly_exp_series(arb_poly_t f, const arb_poly_t h, long n, long prec)
    void _arb_poly_sin_cos_series_basecase(arb_ptr s, arb_ptr c, arb_srcptr h, long hlen, long n, long prec)
    void arb_poly_sin_cos_series_basecase(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)
    void _arb_poly_sin_cos_series_tangent(arb_ptr s, arb_ptr c, const arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_sin_cos_series_tangent(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)
    void _arb_poly_sin_cos_series(arb_ptr s, arb_ptr c, const arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_sin_cos_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)
    void _arb_poly_sin_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)
    void arb_poly_sin_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    void _arb_poly_cos_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)
    void arb_poly_cos_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    void _arb_poly_tan_series(arb_ptr g, arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_tan_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    void _arb_poly_compose_series_brent_kung(arb_ptr res, arb_srcptr poly1, long len1, arb_srcptr poly2, long len2, long n, long prec)
    void arb_poly_compose_series_brent_kung(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, long n, long prec)
    void _arb_poly_evaluate_acb_horner(acb_t res, arb_srcptr f, long len, const acb_t x, long prec)
    void arb_poly_evaluate_acb_horner(acb_t res, const arb_poly_t f, const acb_t a, long prec)
    void _arb_poly_evaluate_acb_rectangular(acb_t y, arb_srcptr poly, long len, const acb_t x, long prec)
    void arb_poly_evaluate_acb_rectangular(acb_t res, const arb_poly_t f, const acb_t a, long prec)
    void _arb_poly_evaluate_acb(acb_t res, arb_srcptr f, long len, const acb_t x, long prec)
    void arb_poly_evaluate_acb(acb_t res, const arb_poly_t f, const acb_t a, long prec)
    void _arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec)
    void arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec)
    void _arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec)
    void arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec)
    void _arb_poly_evaluate2_acb(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec)
    void arb_poly_evaluate2_acb(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec)
    void _arb_poly_gamma_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_gamma_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    void _arb_poly_rgamma_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_rgamma_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    void _arb_poly_lgamma_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_lgamma_series(arb_poly_t res, const arb_poly_t f, long n, long prec)
    void _arb_poly_rising_ui_series(arb_ptr res, arb_srcptr f, long flen, ulong r, long trunc, long prec)
    void arb_poly_rising_ui_series(arb_poly_t res, const arb_poly_t f, ulong r, long trunc, long prec)
    void _arb_poly_zeta_series(arb_ptr res, arb_srcptr h, long hlen, const arb_t a, int deflate, long len, long prec)
    void arb_poly_zeta_series(arb_poly_t res, const arb_poly_t f, const arb_t a, int deflate, long n, long prec)
    void _arb_poly_riemann_siegel_theta_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_riemann_siegel_theta_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void _arb_poly_riemann_siegel_z_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_riemann_siegel_z_series(arb_poly_t res, const arb_poly_t h, long n, long prec)

    void arb_poly_swinnerton_dyer_ui(arb_poly_t poly, ulong n, long prec)
    int arb_poly_get_unique_fmpz_poly(fmpz_poly_t res, const arb_poly_t src)

    void _arb_poly_sin_cos_pi_series(arb_ptr s, arb_ptr c, const arb_srcptr h, long hlen, long len, long prec)
    void arb_poly_sin_cos_pi_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, long n, long prec)
    void _arb_poly_sin_pi_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)
    void arb_poly_sin_pi_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    void _arb_poly_cos_pi_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)
    void arb_poly_cos_pi_series(arb_poly_t g, const arb_poly_t h, long n, long prec)
    void _arb_poly_cot_pi_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec)
    void arb_poly_cot_pi_series(arb_poly_t g, const arb_poly_t h, long n, long prec)

    void arb_poly_lambertw_series(arb_poly_t res, const arb_poly_t z, int flags, long len, long prec)

cdef extern from "arb_mat.h":
    ctypedef struct arb_mat_struct:
        arb_ptr entries
        long r
        long c
        arb_ptr * rows

    ctypedef arb_mat_struct arb_mat_t[1]

    arb_struct * arb_mat_entry(arb_mat_t mat, long i, long j)

    long arb_mat_nrows(const arb_mat_t x)
    long arb_mat_ncols(const arb_mat_t x)

    void arb_mat_init(arb_mat_t mat, long r, long c)
    void arb_mat_clear(arb_mat_t mat)

    void arb_mat_set(arb_mat_t dest, const arb_mat_t src)
    void arb_mat_set_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src)
    void arb_mat_set_fmpq_mat(arb_mat_t dest, const fmpq_mat_t src, long prec)
    void arb_mat_printd(const arb_mat_t mat, long digits)
    int arb_mat_equal(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_overlaps(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_contains(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_contains_fmpq_mat(const arb_mat_t mat1, const fmpq_mat_t mat2)
    int arb_mat_contains_fmpz_mat(const arb_mat_t mat1, const fmpz_mat_t mat2)

    void arb_mat_zero(arb_mat_t mat)
    void arb_mat_one(arb_mat_t mat)

    void arb_mat_bound_inf_norm(mag_t b, const arb_mat_t A)

    void arb_mat_neg(arb_mat_t dest, const arb_mat_t src)
    void arb_mat_add(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    void arb_mat_sub(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    void arb_mat_mul(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec)
    void arb_mat_pow_ui(arb_mat_t B, const arb_mat_t A, ulong exp, long prec)

    void arb_mat_scalar_mul_2exp_si(arb_mat_t B, const arb_mat_t A, long c)
    void arb_mat_scalar_addmul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)
    void arb_mat_scalar_mul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)
    void arb_mat_scalar_div_si(arb_mat_t B, const arb_mat_t A, long c, long prec)
    void arb_mat_scalar_addmul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)
    void arb_mat_scalar_mul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)
    void arb_mat_scalar_div_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)
    void arb_mat_scalar_addmul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
    void arb_mat_scalar_mul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
    void arb_mat_scalar_div_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)

    int arb_mat_lu(long * P, arb_mat_t LU, const arb_mat_t A, long prec)
    void arb_mat_solve_lu_precomp(arb_mat_t X, const long * perm, const arb_mat_t A, const arb_mat_t B, long prec)
    int arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
    int arb_mat_solve_lu(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
    int arb_mat_solve_precond(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)
    int arb_mat_inv(arb_mat_t X, const arb_mat_t A, long prec)
    void arb_mat_det(arb_t det, const arb_mat_t A, long prec)

    void arb_mat_exp(arb_mat_t B, const arb_mat_t A, long prec)

    void _arb_mat_charpoly(arb_ptr cp, const arb_mat_t mat, long prec)
    void arb_mat_charpoly(arb_poly_t cp, const arb_mat_t mat, long prec)

    void arb_mat_transpose(arb_mat_t B, const arb_mat_t A)

    void arb_mat_trace(arb_t trace, const arb_mat_t mat, long prec)
    void arb_mat_ones(arb_mat_t mat)
    void arb_mat_hilbert(arb_mat_t mat, long prec)
    void arb_mat_pascal(arb_mat_t mat, int triangular, long prec)
    void arb_mat_stirling(arb_mat_t mat, int kind, long prec)
    void arb_mat_dct(arb_mat_t mat, int type, long prec)

    void arb_mat_get_mid(arb_mat_t B, const arb_mat_t A)

    int arb_mat_eq(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_ne(const arb_mat_t mat1, const arb_mat_t mat2)
    int arb_mat_equal(const arb_mat_t mat1, const arb_mat_t mat2)

    void arb_mat_frobenius_norm(arb_t res, const arb_mat_t A, long prec)

    int arb_mat_approx_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec)

cdef extern from "acb_poly.h":
    ctypedef struct acb_poly_struct:
        acb_ptr coeffs
        long length
        long alloc

    ctypedef acb_poly_struct acb_poly_t[1]

    void acb_poly_init(acb_poly_t poly)
    void acb_poly_init2(acb_poly_t poly, long len)
    void acb_poly_clear(acb_poly_t poly)
    void acb_poly_fit_length(acb_poly_t poly, long len)
    void _acb_poly_set_length(acb_poly_t poly, long len)
    void _acb_poly_normalise(acb_poly_t poly)
    void acb_poly_swap(acb_poly_t poly1, acb_poly_t poly2)
    long acb_poly_length(const acb_poly_t poly)
    long acb_poly_degree(const acb_poly_t poly)
    void acb_poly_zero(acb_poly_t poly)
    void acb_poly_one(acb_poly_t poly)
    void acb_poly_set_coeff_si(acb_poly_t poly, long n, long x)
    void acb_poly_set_coeff_acb(acb_poly_t poly, long n, const acb_t x)
    void acb_poly_get_coeff_acb(acb_t x, const acb_poly_t poly, long n)
    acb_ptr acb_poly_get_coeff_ptr(arb_poly_t poly, long n)
    void _acb_poly_shift_right(acb_ptr res, acb_srcptr poly, long len, long n)
    void acb_poly_shift_right(acb_poly_t res, const acb_poly_t poly, long n)
    void _acb_poly_shift_left(acb_ptr res, acb_srcptr poly, long len, long n)
    void acb_poly_shift_left(acb_poly_t res, const acb_poly_t poly, long n)
    void acb_poly_truncate(acb_poly_t poly, long newlen)
    void acb_poly_printd(const acb_poly_t poly, long digits)
    void _acb_poly_evaluate_horner(acb_t res, acb_srcptr f, long len, const acb_t a, long prec)
    void acb_poly_evaluate_horner(acb_t res, const acb_poly_t f, const acb_t a, long prec)
    void _acb_poly_evaluate_rectangular(acb_t y, acb_srcptr poly, long len, const acb_t x, long prec)
    void acb_poly_evaluate_rectangular(acb_t res, const acb_poly_t f, const acb_t a, long prec)
    void _acb_poly_evaluate(acb_t res, acb_srcptr f, long len, const acb_t a, long prec)
    void acb_poly_evaluate(acb_t res, const acb_poly_t f, const acb_t a, long prec)
    void _acb_poly_evaluate2_horner(acb_t y, acb_t z, acb_srcptr f, long len, const acb_t x, long prec)
    void acb_poly_evaluate2_horner(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, long prec)
    void _acb_poly_evaluate2_rectangular(acb_t y, acb_t z, acb_srcptr f, long len, const acb_t x, long prec)
    void acb_poly_evaluate2_rectangular(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, long prec)
    void _acb_poly_evaluate2(acb_t y, acb_t z, acb_srcptr f, long len, const acb_t x, long prec)
    void acb_poly_evaluate2(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, long prec)
    void _acb_poly_derivative(acb_ptr res, acb_srcptr poly, long len, long prec)
    void acb_poly_derivative(acb_poly_t res, const acb_poly_t poly, long prec)
    void _acb_poly_integral(acb_ptr res, acb_srcptr poly, long len, long prec)
    void acb_poly_integral(acb_poly_t res, const acb_poly_t poly, long prec)
    void acb_poly_set(acb_poly_t dest, const acb_poly_t src)
    void acb_poly_set_round(acb_poly_t dest, const acb_poly_t src, long prec)
    void acb_poly_set_arb_poly(acb_poly_t poly, const arb_poly_t re)
    void acb_poly_set2_arb_poly(acb_poly_t poly, const arb_poly_t re, const arb_poly_t im)
    void acb_poly_set_fmpq_poly(acb_poly_t poly, const fmpq_poly_t re, long prec)
    void acb_poly_set2_fmpq_poly(acb_poly_t poly, const fmpq_poly_t re, const fmpq_poly_t im, long prec)
    void acb_poly_set_fmpz_poly(acb_poly_t poly, const fmpz_poly_t src, long prec)
    void acb_poly_set_acb(acb_poly_t poly, const acb_t c)
    void acb_poly_set_si(acb_poly_t poly, long c)
    void acb_poly_randtest(acb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits)
    int acb_poly_equal(const acb_poly_t A, const acb_poly_t B)
    int acb_poly_contains_fmpz_poly(const acb_poly_t poly1, const fmpz_poly_t poly2)
    int acb_poly_contains_fmpq_poly(const acb_poly_t poly1, const fmpq_poly_t poly2)
    int _acb_poly_overlaps(acb_srcptr poly1, long len1, acb_srcptr poly2, long len2)
    int acb_poly_overlaps(const acb_poly_t poly1, const acb_poly_t poly2)
    int acb_poly_contains(const acb_poly_t poly1, const acb_poly_t poly2)
    void _acb_poly_add(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long prec)
    void acb_poly_add(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long prec)
    void _acb_poly_sub(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long prec)
    void acb_poly_sub(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long prec)
    void acb_poly_neg(acb_poly_t res, const acb_poly_t poly)
    void acb_poly_scalar_mul_2exp_si(acb_poly_t res, const acb_poly_t poly, long c)
    void acb_poly_mullow_classical(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)
    void _acb_poly_mullow_classical(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)
    void _acb_poly_mullow_transpose(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)
    void acb_poly_mullow_transpose(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)
    void _acb_poly_mullow_transpose_gauss(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)
    void acb_poly_mullow_transpose_gauss(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)
    void _acb_poly_mullow(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)
    void acb_poly_mullow(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)
    void _acb_poly_mul(acb_ptr C, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)
    void acb_poly_mul(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long prec)
    void _acb_poly_inv_series(acb_ptr Qinv, acb_srcptr Q, long Qlen, long len, long prec)
    void acb_poly_inv_series(acb_poly_t Qinv, const acb_poly_t Q, long n, long prec)
    void  _acb_poly_div_series(acb_ptr Q, acb_srcptr A, long Alen, acb_srcptr B, long Blen, long n, long prec)
    void acb_poly_div_series(acb_poly_t Q, const acb_poly_t A, const acb_poly_t B, long n, long prec)
    void _acb_poly_reverse(acb_ptr res, acb_srcptr poly, long len, long n)
    void _acb_poly_div(acb_ptr Q, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)
    void _acb_poly_divrem(acb_ptr Q, acb_ptr R, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)
    void _acb_poly_rem(acb_ptr R, acb_srcptr A, long lenA, acb_srcptr B, long lenB, long prec)
    int acb_poly_divrem(acb_poly_t Q, acb_poly_t R, const acb_poly_t A, const acb_poly_t B, long prec)
    void _acb_poly_div_root(acb_ptr Q, acb_t R, acb_srcptr A, long len, const acb_t c, long prec)
    void _acb_poly_compose(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long prec)
    void acb_poly_compose(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long prec)
    void _acb_poly_compose_horner(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long prec)
    void acb_poly_compose_horner(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long prec)
    void _acb_poly_compose_divconquer(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long prec)
    void acb_poly_compose_divconquer(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long prec)
    void _acb_poly_compose_series_horner(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)
    void acb_poly_compose_series_horner(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)
    void _acb_poly_compose_series_brent_kung(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)
    void acb_poly_compose_series_brent_kung(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)
    void _acb_poly_compose_series(acb_ptr res, acb_srcptr poly1, long len1, acb_srcptr poly2, long len2, long n, long prec)
    void acb_poly_compose_series(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, long n, long prec)
    void _acb_poly_revert_series_lagrange(acb_ptr Qinv, acb_srcptr Q, long Qlen, long n, long prec)
    void acb_poly_revert_series_lagrange(acb_poly_t Qinv, const acb_poly_t Q, long n, long prec)
    void _acb_poly_revert_series_newton(acb_ptr Qinv, acb_srcptr Q, long Qlen, long n, long prec)
    void acb_poly_revert_series_newton(acb_poly_t Qinv, const acb_poly_t Q, long n, long prec)
    void _acb_poly_revert_series_lagrange_fast(acb_ptr Qinv, acb_srcptr Q, long Qlen, long n, long prec)
    void acb_poly_revert_series_lagrange_fast(acb_poly_t Qinv, const acb_poly_t Q, long n, long prec)
    void _acb_poly_revert_series(acb_ptr Qinv, acb_srcptr Q, long Qlen, long n, long prec)
    void acb_poly_revert_series(acb_poly_t Qinv, const acb_poly_t Q, long n, long prec)
    void _acb_poly_evaluate_vec_fast_precomp(acb_ptr vs, acb_srcptr poly, long plen, acb_ptr * tree, long len, long prec)
    void _acb_poly_evaluate_vec_fast(acb_ptr ys, acb_srcptr poly, long plen, acb_srcptr xs, long n, long prec)
    void acb_poly_evaluate_vec_fast(acb_ptr ys, const acb_poly_t poly, acb_srcptr xs, long n, long prec)
    void _acb_poly_evaluate_vec_iter(acb_ptr ys, acb_srcptr poly, long plen, acb_srcptr xs, long n, long prec)
    void acb_poly_evaluate_vec_iter(acb_ptr ys, const acb_poly_t poly, acb_srcptr xs, long n, long prec)
    void _acb_poly_interpolate_barycentric(acb_ptr poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)
    void acb_poly_interpolate_barycentric(acb_poly_t poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)
    void _acb_poly_interpolation_weights(acb_ptr w, acb_ptr * tree, long len, long prec)
    void _acb_poly_interpolate_fast_precomp(acb_ptr poly, acb_srcptr ys, acb_ptr * tree, acb_srcptr weights, long len, long prec)
    void _acb_poly_interpolate_fast(acb_ptr poly, acb_srcptr xs, acb_srcptr ys, long len, long prec)
    void acb_poly_interpolate_fast(acb_poly_t poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)
    void _acb_poly_interpolate_newton(acb_ptr poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)
    void acb_poly_interpolate_newton(acb_poly_t poly, acb_srcptr xs, acb_srcptr ys, long n, long prec)
    void _acb_poly_product_roots(acb_ptr poly, acb_srcptr xs, long n, long prec)
    void acb_poly_product_roots(acb_poly_t poly, acb_srcptr xs, long n, long prec)
    acb_ptr * _acb_poly_tree_alloc(long len)
    void _acb_poly_tree_free(acb_ptr * tree, long len)
    void _acb_poly_tree_build(acb_ptr * tree, acb_srcptr roots, long len, long prec)
    void _acb_poly_root_inclusion(acb_t r, const acb_t m, acb_srcptr poly, acb_srcptr polyder, long len, long prec)
    long _acb_poly_validate_roots(acb_ptr roots, acb_srcptr poly, long len, long prec)
    void _acb_poly_refine_roots_durand_kerner(acb_ptr roots, acb_srcptr poly, long len, long prec)
    long _acb_poly_find_roots(acb_ptr roots, acb_srcptr poly, acb_srcptr initial, long len, long maxiter, long prec)
    long acb_poly_find_roots(acb_ptr roots, const acb_poly_t poly, acb_srcptr initial, long maxiter, long prec)
    void _acb_poly_pow_ui_trunc_binexp(acb_ptr res, acb_srcptr f, long flen, ulong exp, long len, long prec)
    void acb_poly_pow_ui_trunc_binexp(acb_poly_t res, const acb_poly_t poly, ulong exp, long len, long prec)
    void _acb_poly_pow_ui(acb_ptr res, acb_srcptr f, long flen, ulong exp, long prec)
    void acb_poly_pow_ui(acb_poly_t res, const acb_poly_t poly, ulong exp, long prec)
    void _acb_poly_rsqrt_series(acb_ptr g, acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_rsqrt_series(acb_poly_t g, const acb_poly_t h, long n, long prec)
    void _acb_poly_sqrt_series(acb_ptr g, acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_sqrt_series(acb_poly_t g, const acb_poly_t h, long n, long prec)
    void _acb_poly_log_series(acb_ptr res, acb_srcptr f, long flen, long n, long prec)
    void acb_poly_log_series(acb_poly_t res, const acb_poly_t f, long n, long prec)
    void _acb_poly_atan_series(acb_ptr res, acb_srcptr f, long flen, long n, long prec)
    void acb_poly_atan_series(acb_poly_t res, const acb_poly_t f, long n, long prec)
    void _acb_poly_exp_series_basecase(acb_ptr f, acb_srcptr h, long hlen, long n, long prec)
    void acb_poly_exp_series_basecase(acb_poly_t f, const acb_poly_t h, long n, long prec)
    void _acb_poly_exp_series(acb_ptr f, acb_srcptr h, long hlen, long n, long prec)
    void acb_poly_exp_series(acb_poly_t f, const acb_poly_t h, long n, long prec)
    void _acb_poly_sin_cos_series_basecase(acb_ptr s, acb_ptr c, acb_srcptr h, long hlen, long n, long prec)
    void acb_poly_sin_cos_series_basecase(acb_poly_t s, acb_poly_t c, const acb_poly_t h, long n, long prec)
    void _acb_poly_sin_cos_series_tangent(acb_ptr s, acb_ptr c, const acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_sin_cos_series_tangent(acb_poly_t s, acb_poly_t c, const acb_poly_t h, long n, long prec)
    void _acb_poly_sin_cos_series(acb_ptr s, acb_ptr c, const acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_sin_cos_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, long n, long prec)
    void _acb_poly_sin_series(acb_ptr g, acb_srcptr h, long hlen, long n, long prec)
    void acb_poly_sin_series(acb_poly_t g, const acb_poly_t h, long n, long prec)
    void _acb_poly_cos_series(acb_ptr g, acb_srcptr h, long hlen, long n, long prec)
    void acb_poly_cos_series(acb_poly_t g, const acb_poly_t h, long n, long prec)
    void _acb_poly_tan_series(acb_ptr g, acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_tan_series(acb_poly_t g, const acb_poly_t h, long n, long prec)
    void _acb_poly_gamma_series(acb_ptr res, acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_gamma_series(acb_poly_t res, const acb_poly_t f, long n, long prec)
    void _acb_poly_rgamma_series(acb_ptr res, acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_rgamma_series(acb_poly_t res, const acb_poly_t f, long n, long prec)
    void _acb_poly_lgamma_series(acb_ptr res, acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_lgamma_series(acb_poly_t res, const acb_poly_t f, long n, long prec)
    void _acb_poly_rising_ui_series(acb_ptr res, acb_srcptr f, long flen, ulong r, long trunc, long prec)
    void acb_poly_rising_ui_series(acb_poly_t res, const acb_poly_t f, ulong r, long trunc, long prec)
    void _acb_poly_zeta_series(acb_ptr res, acb_srcptr h, long hlen, const acb_t a, int deflate, long len, long prec)
    void acb_poly_zeta_series(acb_poly_t res, const acb_poly_t f, const acb_t a, int deflate, long n, long prec)
    void _acb_poly_polylog_cpx_zeta(acb_ptr w, const acb_t s, const acb_t z, long len, long prec)
    void _acb_poly_polylog_cpx_small(acb_ptr w, const acb_t s, const acb_t z, long len, long prec)
    void _acb_poly_polylog_cpx(acb_ptr w, const acb_t s, const acb_t z, long len, long prec)
    void _acb_poly_polylog_series(acb_ptr res, acb_srcptr s, long slen, const acb_t z, long len, long prec)
    void acb_poly_polylog_series(acb_poly_t res, const acb_poly_t s, const acb_t z, long n, long prec)

    void _acb_poly_pow_series(acb_ptr h, acb_srcptr f, long flen, acb_srcptr g, long glen, long len, long prec)
    void acb_poly_pow_series(acb_poly_t h, const acb_poly_t f, const acb_poly_t g, long len, long prec)
    void _acb_poly_pow_acb_series(acb_ptr h, acb_srcptr f, long flen, const acb_t g, long len, long prec)
    void acb_poly_pow_acb_series(acb_poly_t h, const acb_poly_t f, const acb_t g, long len, long prec)

    void _acb_poly_agm1_series(acb_ptr res, acb_srcptr z, long zlen, long len, long prec)
    void acb_poly_agm1_series(acb_poly_t res, const acb_poly_t z, long n, long prec)

    void _acb_poly_elliptic_k_series(acb_ptr res, acb_srcptr z, long zlen, long len, long prec)
    void acb_poly_elliptic_k_series(acb_poly_t res, const acb_poly_t z, long n, long prec)
    void _acb_poly_elliptic_p_series(acb_ptr res, acb_srcptr z, long zlen, const acb_t tau, long len, long prec)
    void acb_poly_elliptic_p_series(acb_poly_t res, const acb_poly_t z, const acb_t tau, long n, long prec)

    void _acb_poly_erf_series(acb_ptr res, acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_erf_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    int acb_poly_get_unique_fmpz_poly(fmpz_poly_t res, const acb_poly_t src)

    void _acb_poly_sin_cos_pi_series(acb_ptr s, acb_ptr c, const acb_srcptr h, long hlen, long len, long prec)
    void acb_poly_sin_cos_pi_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, long n, long prec)
    void _acb_poly_sin_pi_series(acb_ptr g, acb_srcptr h, long hlen, long n, long prec)
    void acb_poly_sin_pi_series(acb_poly_t g, const acb_poly_t h, long n, long prec)
    void _acb_poly_cos_pi_series(acb_ptr g, acb_srcptr h, long hlen, long n, long prec)
    void acb_poly_cos_pi_series(acb_poly_t g, const acb_poly_t h, long n, long prec)

    void _acb_poly_cot_pi_series(acb_ptr g, acb_srcptr h, long hlen, long n, long prec)
    void acb_poly_cot_pi_series(acb_poly_t g, const acb_poly_t h, long n, long prec)

    void acb_poly_root_bound_fujiwara(mag_t bound, acb_poly_t poly)

    void acb_poly_lambertw_series(acb_poly_t res, const acb_poly_t z, const fmpz_t k, int flags, long len, long prec)

cdef extern from "acb_mat.h":
    ctypedef struct acb_mat_struct:
        acb_ptr entries
        long r
        long c
        acb_ptr * rows

    ctypedef acb_mat_struct acb_mat_t[1]

    acb_struct * acb_mat_entry(acb_mat_t mat, long i, long j)

    long acb_mat_nrows(const acb_mat_t x)
    long acb_mat_ncols(const acb_mat_t x)

    void acb_mat_init(acb_mat_t mat, long r, long c)
    void acb_mat_clear(acb_mat_t mat)

    void acb_mat_set(acb_mat_t dest, const acb_mat_t src)
    void acb_mat_set_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src)
    void acb_mat_set_fmpq_mat(acb_mat_t dest, const fmpq_mat_t src, long prec)
    void acb_mat_printd(const acb_mat_t mat, long digits)
    int acb_mat_equal(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_overlaps(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_contains(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_contains_fmpq_mat(const acb_mat_t mat1, const fmpq_mat_t mat2)
    int acb_mat_contains_fmpz_mat(const acb_mat_t mat1, const fmpz_mat_t mat2)

    void acb_mat_zero(acb_mat_t mat)
    void acb_mat_one(acb_mat_t mat)

    void acb_mat_bound_inf_norm(mag_t b, const acb_mat_t A)

    void acb_mat_neg(acb_mat_t dest, const acb_mat_t src)
    void acb_mat_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_pow_ui(acb_mat_t B, const acb_mat_t A, ulong exp, long prec)

    void acb_mat_scalar_mul_2exp_si(acb_mat_t B, const acb_mat_t A, long c)
    void acb_mat_scalar_addmul_si(acb_mat_t B, const acb_mat_t A, long c, long prec)
    void acb_mat_scalar_mul_si(acb_mat_t B, const acb_mat_t A, long c, long prec)
    void acb_mat_scalar_div_si(acb_mat_t B, const acb_mat_t A, long c, long prec)
    void acb_mat_scalar_addmul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)
    void acb_mat_scalar_mul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)
    void acb_mat_scalar_div_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)
    void acb_mat_scalar_addmul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)
    void acb_mat_scalar_mul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)
    void acb_mat_scalar_div_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)

    int acb_mat_lu(long * P, acb_mat_t LU, const acb_mat_t A, long prec)
    void acb_mat_solve_lu_precomp(acb_mat_t X, const long * perm, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_solve_lu(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_solve_precond(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_inv(acb_mat_t X, const acb_mat_t A, long prec)
    void acb_mat_det(acb_t det, const acb_mat_t A, long prec)

    void acb_mat_exp(acb_mat_t B, const acb_mat_t A, long prec)

    void _acb_mat_charpoly(acb_ptr cp, const acb_mat_t mat, long prec)
    void acb_mat_charpoly(acb_poly_t cp, const acb_mat_t mat, long prec)


    void acb_mat_conjugate(acb_mat_t mat1, const acb_mat_t mat2)
    void acb_mat_transpose(acb_mat_t B, const acb_mat_t A)
    void acb_mat_trace(acb_t trace, const acb_mat_t mat, long prec)
    void acb_mat_get_mid(acb_mat_t B, const acb_mat_t A)

    void acb_mat_dft(acb_mat_t res, int kind, long prec)

    void acb_mat_frobenius_norm(arb_t res, const acb_mat_t A, long prec)

    int acb_mat_eq(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_ne(const acb_mat_t mat1, const acb_mat_t mat2)
    int acb_mat_equal(const acb_mat_t mat1, const acb_mat_t mat2)

    int acb_mat_approx_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)

    void acb_mat_randtest_eig(acb_mat_t A, flint_rand_t state, acb_srcptr E, long prec)
    int acb_mat_approx_eig_qr(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, const mag_t tol, long maxiter, long prec)
    void acb_mat_eig_global_enclosure(mag_t eps, const acb_mat_t A, acb_srcptr E, const acb_mat_t R, long prec)

    int acb_mat_eig_simple_rump(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)
    int acb_mat_eig_simple_vdhoeven_mourrain(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)
    int acb_mat_eig_simple(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)

    int acb_mat_eig_multiple_rump(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)
    int acb_mat_eig_multiple(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, long prec)


cdef extern from "acb_modular.h":
    void acb_modular_theta(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau, long prec)
    void acb_modular_theta_jet(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau, long len, long prec)
    void acb_modular_theta_series(acb_poly_t theta1, acb_poly_t theta2, acb_poly_t theta3, acb_poly_t theta4, const acb_poly_t z, const acb_t tau, long len, long prec)
    void acb_modular_eta(acb_t r, const acb_t tau, long prec)
    void acb_modular_j(acb_t r, const acb_t tau, long prec)
    void acb_modular_lambda(acb_t r, const acb_t tau, long prec)
    void acb_modular_delta(acb_t r, const acb_t tau, long prec)
    void acb_modular_eisenstein(acb_ptr r, const acb_t tau, long len, long prec)
    void acb_modular_elliptic_p(acb_t wp, const acb_t z, const acb_t tau, long prec)
    void acb_modular_elliptic_p_zpx(acb_ptr wp, const acb_t z, const acb_t tau, long len, long prec)
    void acb_modular_elliptic_k(acb_t w, const acb_t m, long prec)
    void acb_modular_elliptic_k_cpx(acb_ptr w, const acb_t m, long len, long prec)
    void acb_modular_elliptic_e(acb_t w, const acb_t m, long prec)
    void acb_modular_hilbert_class_poly(fmpz_poly_t res, long D)

cdef extern from "acb_hypgeom.h":
    void acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, long prec)
    void acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, long prec)
    void acb_hypgeom_bessel_i(acb_t res, const acb_t nu, const acb_t z, long prec)
    void acb_hypgeom_bessel_y(acb_t res, const acb_t nu, const acb_t z, long prec)

    void acb_hypgeom_bessel_k_scaled(acb_t res, const acb_t nu, const acb_t z, long prec)
    void acb_hypgeom_bessel_i_scaled(acb_t res, const acb_t nu, const acb_t z, long prec)

    void acb_hypgeom_erf(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long n, long prec)
    void acb_hypgeom_u_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, long n, long prec)
    void acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t z, long prec)
    void acb_hypgeom_m(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)
    void acb_hypgeom_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)

    long acb_hypgeom_pfq_choose_n(acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long prec)
    void acb_hypgeom_pfq(acb_t res, acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, int regularized, long prec)
    void acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z, int modified, long prec)
    void acb_hypgeom_gamma_upper_asymp(acb_t res, const acb_t s, const acb_t z, int modified, long prec)
    void acb_hypgeom_expint(acb_t res, const acb_t s, const acb_t z, long prec)
    void acb_hypgeom_gamma_lower(acb_t res, const acb_t s, const acb_t z, int modified, long prec)
    void acb_hypgeom_beta_lower(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, long prec)
    void acb_hypgeom_erfc(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_erfi(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_pfq_series_direct(acb_poly_t res, const acb_poly_struct * a, long p, const acb_poly_struct * b, long q, const acb_poly_t z, int regularized, long n, long len, long prec)
    void acb_hypgeom_ei(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_si(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_ci(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_shi(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_chi(acb_t res, const acb_t z, long prec)
    void acb_hypgeom_li(acb_t res, const acb_t z, int offset, long prec)
    void acb_hypgeom_2f1(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, long prec)
    void acb_hypgeom_0f1(acb_t res, const acb_t a, const acb_t z, int regularized, long prec)
    void acb_hypgeom_legendre_p(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, long prec)
    void acb_hypgeom_legendre_q(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, long prec)
    void acb_hypgeom_spherical_y(acb_t res, long n, long m, const acb_t theta, const acb_t phi, long prec)
    void acb_hypgeom_jacobi_p(acb_t res, const acb_t n, const acb_t a, const acb_t b, const acb_t z, long prec)
    void acb_hypgeom_gegenbauer_c(acb_t res, const acb_t n, const acb_t m, const acb_t z, long prec)
    void acb_hypgeom_laguerre_l(acb_t res, const acb_t n, const acb_t m, const acb_t z, long prec)
    void acb_hypgeom_hermite_h(acb_t res, const acb_t n, const acb_t z, long prec)
    void acb_hypgeom_chebyshev_t(acb_t res, const acb_t n, const acb_t z, long prec)
    void acb_hypgeom_chebyshev_u(acb_t res, const acb_t n, const acb_t z, long prec)

    void acb_hypgeom_airy_bound(mag_t ai, mag_t aip, mag_t bi, mag_t bip, const acb_t z)
    void acb_hypgeom_airy_asymp(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, long n, long prec)
    void acb_hypgeom_airy_direct(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, long n, long prec)
    void acb_hypgeom_airy(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, long prec)
    void acb_hypgeom_airy_jet(acb_ptr ai, acb_ptr bi, const acb_t z, long len, long prec)
    void _acb_hypgeom_airy_series(acb_ptr ai, acb_ptr ai_prime, acb_ptr bi, acb_ptr bi_prime, acb_srcptr z, long zlen, long len, long prec)
    void acb_hypgeom_airy_series(acb_poly_t ai, acb_poly_t ai_prime, acb_poly_t bi, acb_poly_t bi_prime, const acb_poly_t z, long len, long prec)

    void acb_hypgeom_coulomb(acb_t F, acb_t G, acb_t Hpos, acb_t Hneg, const acb_t l, const acb_t eta, const acb_t z, long prec)
    void acb_hypgeom_coulomb_series(acb_poly_t F, acb_poly_t G, acb_poly_t Hpos, acb_poly_t Hneg, const acb_t l, const acb_t eta, const acb_poly_t z, long len, long prec)

    void acb_hypgeom_erf_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_erfc_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_erfi_series(acb_poly_t res, const acb_poly_t h, long n, long prec)

    void acb_hypgeom_fresnel(acb_t res1, acb_t res2, const acb_t z, int normalized, long prec)
    void acb_hypgeom_fresnel_series(acb_poly_t res1, acb_poly_t res2, const acb_poly_t h, int normalized, long n, long prec)

    void _acb_hypgeom_gamma_upper_series(acb_ptr g, const acb_t s, acb_srcptr h, long hlen, int regularized, long n, long prec)
    void acb_hypgeom_gamma_upper_series(acb_poly_t g, const acb_t s, const acb_poly_t h, int regularized, long n, long prec)

    void _acb_hypgeom_gamma_lower_series(acb_ptr g, const acb_t s, acb_srcptr h, long hlen, int regularized, long n, long prec)
    void acb_hypgeom_gamma_lower_series(acb_poly_t g, const acb_t s, const acb_poly_t h, int regularized, long n, long prec)

    void _acb_hypgeom_beta_lower_series(acb_ptr g, const acb_t s, const acb_t t, acb_srcptr h, long hlen, int regularized, long n, long prec)
    void acb_hypgeom_beta_lower_series(acb_poly_t g, const acb_t s, const acb_t t, const acb_poly_t h, int regularized, long n, long prec)

    void acb_hypgeom_ei_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_si_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_ci_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_shi_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_chi_series(acb_poly_t res, const acb_poly_t h, long n, long prec)
    void acb_hypgeom_li_series(acb_poly_t res, const acb_poly_t h, int offset, long n, long prec)

cdef extern from "arb_hypgeom.h":
    void arb_hypgeom_pfq(arb_t res, arb_srcptr a, long p, arb_srcptr b, long q, const arb_t z, int regularized, long prec)
    void arb_hypgeom_0f1(arb_t res, const arb_t a, const arb_t z, int regularized, long prec)
    void arb_hypgeom_m(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)
    void arb_hypgeom_1f1(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)
    void arb_hypgeom_u(arb_t res, const arb_t a, const arb_t b, const arb_t z, long prec)
    void arb_hypgeom_2f1(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, long prec)

    void arb_hypgeom_erf(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_erf_series(arb_poly_t g, const arb_poly_t h, long len, long prec)
    void arb_hypgeom_erfc(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_erfc_series(arb_poly_t g, const arb_poly_t h, long len, long prec)
    void arb_hypgeom_erfi(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_erfi_series(arb_poly_t g, const arb_poly_t h, long len, long prec)
    void arb_hypgeom_fresnel(arb_t res1, arb_t res2, const arb_t z, int normalized, long prec)
    void arb_hypgeom_fresnel_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, int normalized, long len, long prec)

    void arb_hypgeom_ei(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_si(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_ci(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_shi(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_chi(arb_t res, const arb_t z, long prec)
    void arb_hypgeom_li(arb_t res, const arb_t z, int offset, long prec)
    void arb_hypgeom_ei_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_si_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_ci_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_shi_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_chi_series(arb_poly_t res, const arb_poly_t h, long n, long prec)
    void arb_hypgeom_li_series(arb_poly_t res, const arb_poly_t h, int offset, long n, long prec)

    void arb_hypgeom_bessel_j(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_bessel_k(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_bessel_i(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_bessel_y(arb_t res, const arb_t nu, const arb_t z, long prec)

    void arb_hypgeom_bessel_k_scaled(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_bessel_i_scaled(arb_t res, const arb_t nu, const arb_t z, long prec)

    void arb_hypgeom_airy(arb_t ai, arb_t aip, arb_t bi, arb_t bip, const arb_t z, long prec)
    void arb_hypgeom_airy_series(arb_poly_t ai, arb_poly_t ai_prime, arb_poly_t bi, arb_poly_t bi_prime, const arb_poly_t z, long len, long prec)
    void arb_hypgeom_airy_zero(arb_t ai, arb_t aip, arb_t bi, arb_t bip, const fmpz_t n, long prec)

    void arb_hypgeom_coulomb(arb_t F, arb_t G, const arb_t l, const arb_t eta, const arb_t z, long prec)
    void arb_hypgeom_coulomb_series(arb_poly_t F, arb_poly_t G, const arb_t l, const arb_t eta, const arb_poly_t z, long len, long prec)

    void arb_hypgeom_expint(arb_t res, const arb_t s, const arb_t z, long prec)
    void arb_hypgeom_gamma_upper(arb_t res, const arb_t s, const arb_t z, int modified, long prec)
    void arb_hypgeom_gamma_lower(arb_t res, const arb_t s, const arb_t z, int modified, long prec)
    void arb_hypgeom_beta_lower(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, long prec)

    void arb_hypgeom_gamma_upper_series(arb_poly_t g, const arb_t s, const arb_poly_t h, int regularized, long n, long prec)
    void arb_hypgeom_gamma_lower_series(arb_poly_t g, const arb_t s, const arb_poly_t h, int regularized, long n, long prec)
    void arb_hypgeom_beta_lower_series(arb_poly_t g, const arb_t s, const arb_t t, const arb_poly_t h, int regularized, long n, long prec)

    void arb_hypgeom_chebyshev_t(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_chebyshev_u(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_jacobi_p(arb_t res, const arb_t n, const arb_t a, const arb_t b, const arb_t z, long prec)
    void arb_hypgeom_gegenbauer_c(arb_t res, const arb_t n, const arb_t m, const arb_t z, long prec)
    void arb_hypgeom_laguerre_l(arb_t res, const arb_t n, const arb_t m, const arb_t z, long prec)
    void arb_hypgeom_hermite_h(arb_t res, const arb_t nu, const arb_t z, long prec)
    void arb_hypgeom_legendre_p(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, long prec)
    void arb_hypgeom_legendre_q(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, long prec)

    void arb_hypgeom_legendre_p_ui_root(arb_t res, arb_t weight, ulong n, ulong k, long prec)

cdef extern from "dirichlet.h":
    ctypedef struct dirichlet_group_struct:
        ulong q
        ulong q_even
        nmod_t mod
        ulong rad_q
        ulong phi_q
        long neven
        long num
        ulong expo
        void * P
        ulong * generators
        ulong * PHI
    ctypedef dirichlet_group_struct dirichlet_group_t[1]

    ctypedef struct dirichlet_char_struct:
        ulong n
        ulong * log
    ctypedef dirichlet_char_struct dirichlet_char_t[1]

    ulong dirichlet_group_size(const dirichlet_group_t G)
    void dirichlet_group_init(dirichlet_group_t G, ulong q)
    void dirichlet_group_clear(dirichlet_group_t G)
    ulong dirichlet_number_primitive(const dirichlet_group_t G)

    void dirichlet_char_init(dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_clear(dirichlet_char_t x)
    void dirichlet_char_print(const dirichlet_group_t G, const dirichlet_char_t x)

    void dirichlet_char_set(dirichlet_char_t x, const dirichlet_group_t G, const dirichlet_char_t y)
    int dirichlet_char_eq(const dirichlet_char_t x, const dirichlet_char_t y)
    int dirichlet_parity_char(const dirichlet_group_t G, const dirichlet_char_t x)
    ulong dirichlet_conductor_char(const dirichlet_group_t G, const dirichlet_char_t x)
    ulong dirichlet_order_char(const dirichlet_group_t G, const dirichlet_char_t x)

    void dirichlet_char_log(dirichlet_char_t x, const dirichlet_group_t G, ulong m)
    ulong dirichlet_char_exp(const dirichlet_group_t G, const dirichlet_char_t x)
    ulong _dirichlet_char_exp(dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_index(dirichlet_char_t x, const dirichlet_group_t G, ulong j)
    ulong dirichlet_index_char(const dirichlet_group_t G, const dirichlet_char_t x)
    void dirichlet_char_one(dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_first_primitive(dirichlet_char_t x, const dirichlet_group_t G)
    int dirichlet_char_next(dirichlet_char_t x, const dirichlet_group_t G)
    int dirichlet_char_next_primitive(dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_mul(dirichlet_char_t c, const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b)
    void dirichlet_char_pow(dirichlet_char_t c, const dirichlet_group_t G, const dirichlet_char_t a, ulong n)
    void dirichlet_char_lower(dirichlet_char_t y, const dirichlet_group_t H, const dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_lift(dirichlet_char_t y, const dirichlet_group_t G, const dirichlet_char_t x, const dirichlet_group_t H)

    cdef ulong DIRICHLET_CHI_NULL

    ulong dirichlet_pairing(const dirichlet_group_t G, ulong m, ulong n)
    ulong dirichlet_pairing_char(const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b)

    int dirichlet_char_is_principal(const dirichlet_group_t G, const dirichlet_char_t chi)
    int dirichlet_char_is_real(const dirichlet_group_t G, const dirichlet_char_t chi)
    int dirichlet_char_is_primitive(const dirichlet_group_t G, const dirichlet_char_t chi)
    ulong dirichlet_chi(const dirichlet_group_t G, const dirichlet_char_t chi, ulong n)

cdef extern from "acb_dirichlet.h":
    void acb_dirichlet_eta(acb_t res, const acb_t s, long prec)
    void acb_dirichlet_chi(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, ulong n, long prec)

    void acb_dirichlet_l(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, long prec)
    void acb_dirichlet_hardy_z(acb_ptr res, const acb_t t, const dirichlet_group_t G, const dirichlet_char_t chi, long len, long prec)
    void acb_dirichlet_l_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, long len, long prec)

    void acb_dirichlet_stieltjes(acb_t res, const fmpz_t n, const acb_t a, long prec)

    void acb_dirichlet_gram_point(arb_t res, const fmpz_t n, const dirichlet_group_t G, const dirichlet_char_t chi, long prec)
    void acb_dirichlet_zeta_zeros(acb_ptr res, const fmpz_t n, long len, long prec)
    void acb_dirichlet_zeta_nzeros(arb_t res, const arb_t t, long prec)
    void acb_dirichlet_backlund_s(arb_t res, const arb_t t, long prec)
    void acb_dirichlet_zeta_zero(acb_t res, const fmpz_t n, long prec)
    void acb_dirichlet_zeta_zeros(acb_ptr res, const fmpz_t n, long len, long prec)

cdef extern from "acb_elliptic.h":
    void acb_elliptic_rf(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, long prec)
    void acb_elliptic_rj(acb_t res, const acb_t x, const acb_t y, const acb_t z, const acb_t p, int flags, long prec)
    void acb_elliptic_rg(acb_t res, const acb_t x, const acb_t y, const acb_t z, int flags, long prec)
    void acb_elliptic_f(acb_t res, const acb_t phi, const acb_t m, int times_pi, long prec)
    void acb_elliptic_e_inc(acb_t res, const acb_t phi, const acb_t m, int times_pi, long prec)
    void acb_elliptic_pi(acb_t res, const acb_t n, const acb_t m, long prec)
    void acb_elliptic_pi_inc(acb_t res, const acb_t n, const acb_t phi, const acb_t m, int times_pi, long prec)
    void acb_elliptic_p(acb_t res, const acb_t z, const acb_t tau, long prec)
    void acb_elliptic_zeta(acb_t res, const acb_t z, const acb_t tau, long prec)
    void acb_elliptic_sigma(acb_t res, const acb_t z, const acb_t tau, long prec)
    void acb_elliptic_roots(acb_t e1, acb_t e2, acb_t e3, const acb_t tau, long prec)
    void acb_elliptic_invariants(acb_t g2, acb_t g3, const acb_t tau, long prec)
    void acb_elliptic_inv_p(acb_t res, const acb_t z, const acb_t tau, long prec)

cdef extern from "acb_calc.h":
    ctypedef int (*acb_calc_func_t)(acb_ptr out, const acb_t inp, void * param, long order, long prec)

    ctypedef struct acb_calc_integrate_opt_struct:
        long deg_limit
        long eval_limit
        long depth_limit
        int use_heap
        int verbose

    ctypedef acb_calc_integrate_opt_struct acb_calc_integrate_opt_t[1]

    void acb_calc_integrate_opt_init(acb_calc_integrate_opt_t options)

    int acb_calc_integrate(acb_t res, acb_calc_func_t f, void * param,
        const acb_t a, const acb_t b,
        long goal, const mag_t tol,
        const acb_calc_integrate_opt_t options,
        long prec)

cdef extern from "arb_fmpz_poly.h":
    void arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t poly, const arb_t x, long prec)
    void arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t poly, const acb_t x, long prec)
    void arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, long prec)
    ulong arb_fmpz_poly_deflation(const fmpz_poly_t poly)
    void arb_fmpz_poly_deflate(fmpz_poly_t res, const fmpz_poly_t poly, ulong deflation)

cdef extern from "acb_dft.h":
    void acb_dft(acb_ptr w, acb_srcptr v, long n, long prec)
    void acb_dft_inverse(acb_ptr w, acb_srcptr v, long n, long prec)

cdef extern from "flint/mpoly.h":
    ctypedef enum ordering_t:
        ORD_LEX
        ORD_DEGLEX
        ORD_DEGREVLEX

    ctypedef struct mpoly_ctx_struct:
        slong nvars
        slong nfields
        ordering_t ord
        int deg
        int rev
        slong lut_words_per_exp[FLINT_BITS]
        unsigned char lut_fix_bits[FLINT_BITS]

    ctypedef mpoly_ctx_struct mpoly_ctx_t[1]

    void mpoly_ctx_init(mpoly_ctx_t ctx, slong nvars, const ordering_t ord)
    void mpoly_ctx_clear(mpoly_ctx_t mctx)


cdef extern from "flint/fmpz_mpoly.h":
    ctypedef struct fmpz_mpoly_ctx_struct:
        mpoly_ctx_t minfo

    ctypedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1]

    ctypedef struct fmpz_mpoly_struct:
        fmpz_struct * coeffs
        ulong * exps
        slong alloc
        slong length
        flint_bitcnt_t bits

    ctypedef fmpz_mpoly_struct fmpz_mpoly_t[1]

    void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, slong nvars, const ordering_t ord)
    void fmpz_mpoly_ctx_init_rand(fmpz_mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars)
    void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx)
    slong fmpz_mpoly_ctx_nvars(const fmpz_mpoly_ctx_t ctx)
    ordering_t fmpz_mpoly_ctx_ord(const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_init(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_clear(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    int fmpz_mpoly_set_str_pretty(fmpz_mpoly_t A, const char * str, const char ** x, const fmpz_mpoly_ctx_t ctx)
    char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_is_gen(const fmpz_mpoly_t poly, slong k, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_equal(const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_swap(fmpz_mpoly_t A, fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    slong fmpz_mpoly_max_bits(const fmpz_mpoly_t A)

    int fmpz_mpoly_is_fmpz(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_get_fmpz(fmpz_t c, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_fmpz(fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_ui(fmpz_mpoly_t A, ulong c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_si(fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_zero(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_one(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_equal_fmpz(const fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_equal_ui(const fmpz_mpoly_t A, ulong c, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_equal_si(const fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_is_zero(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_is_one(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    int fmpz_mpoly_degrees_fit_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_degrees_fmpz(fmpz_struct ** degs, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_degrees_si(slong * degs, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)
    slong fmpz_mpoly_degree_si(const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_total_degree_fits_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_total_degree_fmpz(fmpz_t td, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    slong fmpz_mpoly_total_degree_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mpoly_t A, const fmpz_mpoly_t M, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_coeff_fmpz_monomial(fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_t M, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t A, fmpz_struct * const * exp, const fmpz_mpoly_ctx_t ctx)
    ulong fmpz_mpoly_get_coeff_ui_fmpz(const fmpz_mpoly_t A, fmpz_struct * const * exp, const fmpz_mpoly_ctx_t ctx)
    slong fmpz_mpoly_get_coeff_si_fmpz(const fmpz_mpoly_t A, fmpz_struct * const * exp, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
    ulong fmpz_mpoly_get_coeff_ui_ui(const fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
    slong fmpz_mpoly_get_coeff_si_ui(const fmpz_mpoly_t A,
                                const ulong * exp, const fmpz_mpoly_ctx_t ctx)
    void _fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t A,
                 const fmpz_t c, const fmpz_struct * exp, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t A,
               const fmpz_t c, fmpz_struct * const * exp, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_coeff_ui_fmpz(fmpz_mpoly_t A,
                const ulong c, fmpz_struct * const * exp, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_coeff_si_fmpz(fmpz_mpoly_t A,
                const slong c, fmpz_struct * const * exp, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_coeff_fmpz_ui(fmpz_mpoly_t A,
                const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_coeff_ui_ui(fmpz_mpoly_t A,
                 const ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_set_coeff_si_ui(fmpz_mpoly_t A,
                 const slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_get_coeff_vars_ui(fmpz_mpoly_t C,
             const fmpz_mpoly_t A,  slong * vars, ulong * exps, slong length,
                                                   const fmpz_mpoly_ctx_t ctx)

    int fmpz_mpoly_cmp(const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

    slong fmpz_mpoly_length(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_get_term_exp_fmpz(fmpz_struct ** exp, const fmpz_mpoly_t A,
                                              slong i, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mpoly_t A,
                                              slong i, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_get_term_exp_si(slong * exp, const fmpz_mpoly_t A,
                                              slong i, const fmpz_mpoly_ctx_t ctx)

    ulong fmpz_mpoly_get_term_var_exp_ui(const fmpz_mpoly_t A, slong i,
                                            slong var, const fmpz_mpoly_ctx_t ctx)

    slong fmpz_mpoly_get_term_var_exp_si(const fmpz_mpoly_t A, slong i,
                                            slong var, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_set_term_exp_fmpz(fmpz_mpoly_t A,
                          slong i, fmpz_struct * const * exp, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_set_term_exp_ui(fmpz_mpoly_t A,
                           slong i, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_get_term(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                              slong i, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_get_term_monomial(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                              slong i, const fmpz_mpoly_ctx_t ctx)

    # Addition/Subtraction
    void fmpz_mpoly_add_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_add_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_sub_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_sub_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_add(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_sub(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

    # Scalar operations
    void fmpz_mpoly_neg(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_scalar_mul_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_scalar_divexact_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_scalar_divexact_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_scalar_divides_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_scalar_divides_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_scalar_divides_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)

    # Differentiation/Integration
    void fmpz_mpoly_derivative(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_integral(fmpz_mpoly_t A, fmpz_t scale, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)

    # Evaluation
    int fmpz_mpoly_evaluate_all_fmpz(fmpz_t ev, const fmpz_mpoly_t A, fmpz_struct * const * vals, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_t val, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B, fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctxB)
    int fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
    void fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_t A, const fmpz_mpoly_t B, const slong * c, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)

    # Multiplication
    void fmpz_mpoly_mul(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

    # Powering
    int fmpz_mpoly_pow_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t k, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_pow_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong k, const fmpz_mpoly_ctx_t ctx)

    # Division
    int fmpz_mpoly_divides(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_divrem(fmpz_mpoly_t Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_quasidivrem(fmpz_t scale, fmpz_mpoly_t Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_div(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_quasidiv(fmpz_t scale, fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B, slong len, const fmpz_mpoly_ctx_t ctx)
    void fmpz_mpoly_term_content(fmpz_mpoly_t M, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
    int fmpz_mpoly_gcd(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

"""
cdef extern from "flint/fmpz_mpoly_factor.h":

    ctypedef struct fmpz_mpoly_factor_struct:
        fmpz_t content
        fmpz_mpoly_struct * poly
        fmpz_struct * exp
        slong length
        slong alloc

    ctypedef fmpz_mpoly_factor_struct fmpz_mpoly_factor_t[1]


    void fmpz_mpoly_factor_init(fmpz_mpoly_factor_t fac, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_t fac, const fmpz_mpoly_ctx_t ctx)

    int fmpz_mpoly_factor(fmpz_mpoly_factor_t fac, const fmpz_mpoly_t A, int full, const fmpz_mpoly_ctx_t ctx)
"""

