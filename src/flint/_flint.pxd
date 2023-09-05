# _flint.pxd
#
# Define the contents of the Python, GMP, Flint and Arb headers.

cdef extern from "Python.h":
    ctypedef void PyObject
#    ctypedef void PyTypeObject
#     ctypedef long Py_ssize_t
#     int PyObject_TypeCheck(object, PyTypeObject*)
#     int PyInt_Check(PyObject *o)
#     int PyLong_Check(PyObject *o)
#     long PyInt_AS_LONG(PyObject *io)
#     double PyFloat_AS_DOUBLE(PyObject *io)
#     Py_ssize_t PyList_GET_SIZE(PyObject *list)
#     long PyLong_AsLongAndOverflow(PyObject *pylong, int *overflow)
#     long long PyLong_AsLongLongAndOverflow(PyObject *pylong, int *overflow)

cdef enum:
    FMPZ_UNKNOWN = 0
    FMPZ_REF = 1
    FMPZ_TMP = 2


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

from flint.flintlib.nmod_vec cimport nmod_t
from flint.flintlib.nmod_poly cimport nmod_poly_t
from flint.flintlib.nmod_mat cimport nmod_mat_t
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.fmpz_poly cimport fmpz_poly_t, fmpz_poly_struct
from flint.flintlib.fmpz_mat cimport fmpz_mat_struct, fmpz_mat_t
from flint.flintlib.fmpq cimport fmpq_t, fmpq_struct
from flint.flintlib.fmpq_poly cimport fmpq_poly_struct, fmpq_poly_t
from flint.flintlib.fmpq_mat cimport fmpq_mat_t
from flint.flintlib.mag cimport mag_struct, mag_t, mag_ptr, mag_srcptr

from flint.flintlib.arf cimport arf_struct, arf_t, arf_ptr, arf_srcptr, arf_rnd_t
from flint.flintlib.arb cimport arb_struct, arb_ptr, arb_srcptr, arb_t
from flint.flintlib.arb cimport arb_midref, arb_radref
from flint.flintlib.acb cimport acb_struct, acb_ptr, acb_srcptr, acb_t
from flint.flintlib.acb cimport acb_realref, acb_imagref
from flint.flintlib.arb_poly cimport arb_poly_struct, arb_poly_t
from flint.flintlib.arb_mat cimport arb_mat_struct, arb_mat_t
from flint.flintlib.acb_poly cimport acb_poly_struct, acb_poly_t
from flint.flintlib.acb_mat cimport acb_mat_struct, acb_mat_t


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

