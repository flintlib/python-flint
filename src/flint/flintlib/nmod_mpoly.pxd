from flint.flintlib.flint cimport fmpz_struct, flint_rand_t, mp_limb_t, slong, ulong, flint_bitcnt_t, nmod_t
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.mpoly cimport mpoly_ctx_t, ordering_t
from flint.flintlib.nmod_poly cimport nmod_poly_struct, nmod_poly_t


cdef extern from "flint/nmod_mpoly.h":
    ctypedef struct nmod_mpoly_struct:
        mp_limb_t * coeffs
        ulong * exps
        slong length
        flint_bitcnt_t bits
        slong coeffs_alloc
        slong exps_alloc

    ctypedef nmod_mpoly_struct nmod_mpoly_t[1]

    ctypedef struct nmod_mpoly_ctx_struct:
        mpoly_ctx_t minfo
        nmod_t mod

    ctypedef nmod_mpoly_ctx_struct nmod_mpoly_ctx_t[1]

    ctypedef struct nmod_mpoly_univar_struct:
        nmod_mpoly_struct * coeffs
        fmpz_struct * exps
        slong alloc
        slong length

    ctypedef nmod_mpoly_univar_struct nmod_mpoly_univar_t[1]

# from here on is parsed
    void nmod_mpoly_ctx_init(nmod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, mp_limb_t n)
    slong nmod_mpoly_ctx_nvars(const nmod_mpoly_ctx_t ctx)
    ordering_t nmod_mpoly_ctx_ord(const nmod_mpoly_ctx_t ctx)
    mp_limb_t nmod_mpoly_ctx_modulus(const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_ctx_clear(nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_init(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_init2(nmod_mpoly_t A, slong alloc, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_init3(nmod_mpoly_t A, slong alloc, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_fit_length(nmod_mpoly_t A, slong len, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_realloc(nmod_mpoly_t A, slong alloc, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_clear(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    char * nmod_mpoly_get_str_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)
    # int nmod_mpoly_fprint_pretty(FILE * file, const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_print_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_set_str_pretty(nmod_mpoly_t A, const char * str, const char ** x, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_gen(nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_is_gen(const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_equal(const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_swap(nmod_mpoly_t A, nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_is_ui(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    ulong nmod_mpoly_get_ui(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set_ui(nmod_mpoly_t A, ulong c, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_zero(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_one(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_equal_ui(const nmod_mpoly_t A, ulong c, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_is_zero(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_is_one(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_degrees_fit_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_degrees_fmpz(fmpz_struct ** degs, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_degrees_si(slong * degs, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_degree_fmpz(fmpz_t deg, const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_degree_si(const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_total_degree_fits_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_total_degree_fmpz(fmpz_t tdeg, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_total_degree_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_used_vars(int * used, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    ulong nmod_mpoly_get_coeff_ui_monomial(const nmod_mpoly_t A, const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set_coeff_ui_monomial(nmod_mpoly_t A, ulong c, const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx)
    ulong nmod_mpoly_get_coeff_ui_fmpz(const nmod_mpoly_t A, fmpz_struct * const * exp, const nmod_mpoly_ctx_t ctx)
    ulong nmod_mpoly_get_coeff_ui_ui(const nmod_mpoly_t A, const ulong * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set_coeff_ui_fmpz(nmod_mpoly_t A, ulong c, fmpz_struct * const * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set_coeff_ui_ui(nmod_mpoly_t A, ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_get_coeff_vars_ui(nmod_mpoly_t C, const nmod_mpoly_t A, const slong * vars, const ulong * exps, slong length, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_cmp(const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    mp_limb_t * nmod_mpoly_term_coeff_ref(nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_is_canonical(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_length(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_resize(nmod_mpoly_t A, slong new_length, const nmod_mpoly_ctx_t ctx)
    ulong nmod_mpoly_get_term_coeff_ui(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set_term_coeff_ui(nmod_mpoly_t A, slong i, ulong c, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_term_exp_fits_si(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_term_exp_fits_ui(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_get_term_exp_fmpz(fmpz_struct ** exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_get_term_exp_ui(ulong * exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_get_term_exp_si(slong * exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    ulong nmod_mpoly_get_term_var_exp_ui(const nmod_mpoly_t A, slong i, slong var, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_get_term_var_exp_si(const nmod_mpoly_t A, slong i, slong var, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set_term_exp_fmpz(nmod_mpoly_t A, slong i, fmpz_struct * const * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_set_term_exp_ui(nmod_mpoly_t A, slong i, const ulong * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_get_term(nmod_mpoly_t M, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_get_term_monomial(nmod_mpoly_t M, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_push_term_ui_fmpz(nmod_mpoly_t A, ulong c, fmpz_struct * const * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_push_term_ui_ffmpz(nmod_mpoly_t A, ulong c, const fmpz_struct * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_push_term_ui_ui(nmod_mpoly_t A, ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_sort_terms(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_combine_like_terms(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_reverse(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_randtest_bound(nmod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_randtest_bounds(nmod_mpoly_t A, flint_rand_t state, slong length, ulong * exp_bounds, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_randtest_bits(nmod_mpoly_t A, flint_rand_t state, slong length, mp_limb_t exp_bits, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_add_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_sub_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_add(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_sub(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_neg(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_make_monic(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_derivative(nmod_mpoly_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
    ulong nmod_mpoly_evaluate_all_ui(const nmod_mpoly_t A, const ulong * vals, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_evaluate_one_ui(nmod_mpoly_t A, const nmod_mpoly_t B, slong var, ulong val, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_compose_nmod_poly(nmod_poly_t A, const nmod_mpoly_t B, nmod_poly_struct * const * C, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_compose_nmod_mpoly_geobucket(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
    int nmod_mpoly_compose_nmod_mpoly_horner(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
    int nmod_mpoly_compose_nmod_mpoly(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
    void nmod_mpoly_compose_nmod_mpoly_gen(nmod_mpoly_t A, const nmod_mpoly_t B, const slong * c, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
    void nmod_mpoly_mul(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_mul_johnson(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_mul_array(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_mul_array_threaded(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_mul_dense(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_pow_fmpz(nmod_mpoly_t A, const nmod_mpoly_t B, const fmpz_t k, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_pow_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong k, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_divides(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_div(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_divrem(nmod_mpoly_t Q, nmod_mpoly_t R, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_divrem_ideal(nmod_mpoly_struct ** Q, nmod_mpoly_t R, const nmod_mpoly_t A, nmod_mpoly_struct * const * B, slong len, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_divides_dense(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_divides_heap_threaded(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_term_content(nmod_mpoly_t M, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_content_vars(nmod_mpoly_t g, const nmod_mpoly_t A, slong * vars, slong vars_length, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_gcd(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_gcd_cofactors(nmod_mpoly_t G, nmod_mpoly_t Abar, nmod_mpoly_t Bbar, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_gcd_brown(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_gcd_hensel(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_resultant(nmod_mpoly_t R, const nmod_mpoly_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_discriminant(nmod_mpoly_t D, const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_sqrt(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_is_square(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_quadratic_root(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_univar_init(nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_univar_clear(nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_univar_swap(nmod_mpoly_univar_t A, nmod_mpoly_univar_t B, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_to_univar(nmod_mpoly_univar_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_from_univar(nmod_mpoly_t A, const nmod_mpoly_univar_t B, slong var, const nmod_mpoly_ctx_t ctx)
    int nmod_mpoly_univar_degree_fits_si(const nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_univar_length(const nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
    slong nmod_mpoly_univar_get_term_exp_si(nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_univar_get_term_coeff(nmod_mpoly_t c, const nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_univar_swap_term_coeff(nmod_mpoly_t c, nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_pow_rmul(nmod_mpoly_t A, const nmod_mpoly_t B, ulong k, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_div_monagan_pearce(nmod_mpoly_t polyq, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_divrem_monagan_pearce(nmod_mpoly_t q, nmod_mpoly_t r, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_divrem_ideal_monagan_pearce(nmod_mpoly_struct ** q, nmod_mpoly_t r, const nmod_mpoly_t poly2, nmod_mpoly_struct * const * poly3, slong len, const nmod_mpoly_ctx_t ctx)
