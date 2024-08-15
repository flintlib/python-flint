from flint.flintlib.flint cimport flint_rand_t, slong, ulong
from flint.flintlib.fq_zech cimport fq_zech_ctx_t, fq_zech_t, fq_zech_struct
from flint.flintlib.fq_zech cimport fq_zech_poly_struct, fq_zech_poly_t

cdef extern from "flint/fq_zech_poly_factor.h":
    # Type definitions **********************************************/
    ctypedef struct fq_zech_poly_factor_struct:
        fq_zech_poly_struct * poly
        slong * exp
        slong num
        slong alloc
    ctypedef fq_zech_poly_factor_struct fq_zech_poly_factor_t[1]

    # Parsed from here **********************************************/
    void fq_zech_poly_factor_init(fq_zech_poly_factor_t fac, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_clear(fq_zech_poly_factor_t fac, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_realloc(fq_zech_poly_factor_t fac, slong alloc, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_fit_length(fq_zech_poly_factor_t fac, slong len, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_set(fq_zech_poly_factor_t res, const fq_zech_poly_factor_t fac, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_print_pretty(const fq_zech_poly_factor_t fac, const char * var, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_print(const fq_zech_poly_factor_t fac, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_insert(fq_zech_poly_factor_t fac, const fq_zech_poly_t poly, slong exp, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_concat(fq_zech_poly_factor_t res, const fq_zech_poly_factor_t fac, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_pow(fq_zech_poly_factor_t fac, slong exp, const fq_zech_ctx_t ctx)
    ulong fq_zech_poly_remove(fq_zech_poly_t f, const fq_zech_poly_t p, const fq_zech_ctx_t ctx)
    int fq_zech_poly_is_irreducible(const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    int fq_zech_poly_is_irreducible_ddf(const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    int fq_zech_poly_is_irreducible_ben_or(const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    int _fq_zech_poly_is_squarefree(const fq_zech_struct * f, slong len, const fq_zech_ctx_t ctx)
    int fq_zech_poly_is_squarefree(const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    int fq_zech_poly_factor_equal_deg_prob(fq_zech_poly_t factor, flint_rand_t state, const fq_zech_poly_t pol, slong d, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_equal_deg(fq_zech_poly_factor_t factors, const fq_zech_poly_t pol, slong d, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_split_single(fq_zech_poly_t linfactor, const fq_zech_poly_t input, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_distinct_deg(fq_zech_poly_factor_t res, const fq_zech_poly_t poly, slong * const * degs, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_squarefree(fq_zech_poly_factor_t res, const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor(fq_zech_poly_factor_t res, fq_zech_t lead, const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_cantor_zassenhaus(fq_zech_poly_factor_t res, const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_kaltofen_shoup(fq_zech_poly_factor_t res, const fq_zech_poly_t poly, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_berlekamp(fq_zech_poly_factor_t factors, const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_with_berlekamp(fq_zech_poly_factor_t res, fq_zech_t leading_coeff, const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_with_cantor_zassenhaus(fq_zech_poly_factor_t res, fq_zech_t leading_coeff, const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    void fq_zech_poly_factor_with_kaltofen_shoup(fq_zech_poly_factor_t res, fq_zech_t leading_coeff, const fq_zech_poly_t f, const fq_zech_ctx_t ctx)
    void fq_zech_poly_iterated_frobenius_preinv(fq_zech_poly_t * rop, slong n, const fq_zech_poly_t v, const fq_zech_poly_t vinv, const fq_zech_ctx_t ctx)
    void fq_zech_poly_roots(fq_zech_poly_factor_t r, const fq_zech_poly_t f, int with_multiplicity, const fq_zech_ctx_t ctx)