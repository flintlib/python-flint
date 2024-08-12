from flint.flintlib.flint cimport slong, flint_rand_t
from flint.flintlib.fq cimport fq_t, fq_struct, fq_ctx_t
from flint.flintlib.fq_poly cimport fq_poly_t



cdef extern from "flint/fq_poly_factor.h":
    # Type definitions **********************************************/
    ctypedef struct fq_poly_factor_struct:
        fq_poly_struct * poly
        slong * exp
        slong num
        slong alloc
    ctypedef fq_poly_factor_struct fq_poly_factor_t[1]

    # Parsed from here **********************************************/
    void fq_poly_factor_init(fq_poly_factor_t fac, const fq_ctx_t ctx)
    void fq_poly_factor_clear(fq_poly_factor_t fac, const fq_ctx_t ctx)
    void fq_poly_factor_realloc(fq_poly_factor_t fac, slong alloc, const fq_ctx_t ctx)
    void fq_poly_factor_fit_length(fq_poly_factor_t fac, slong len, const fq_ctx_t ctx)
    void fq_poly_factor_set(fq_poly_factor_t res, const fq_poly_factor_t fac, const fq_ctx_t ctx)
    void fq_poly_factor_print_pretty(const fq_poly_factor_t fac, const char * var, const fq_ctx_t ctx)
    void fq_poly_factor_print(const fq_poly_factor_t fac, const fq_ctx_t ctx)
    void fq_poly_factor_insert(fq_poly_factor_t fac, const fq_poly_t poly, slong exp, const fq_ctx_t ctx)
    void fq_poly_factor_concat(fq_poly_factor_t res, const fq_poly_factor_t fac, const fq_ctx_t ctx)
    void fq_poly_factor_pow(fq_poly_factor_t fac, slong exp, const fq_ctx_t ctx)
    ulong fq_poly_remove(fq_poly_t f, const fq_poly_t p, const fq_ctx_t ctx)
    int fq_poly_is_irreducible(const fq_poly_t f, const fq_ctx_t ctx)
    int fq_poly_is_irreducible_ddf(const fq_poly_t f, const fq_ctx_t ctx)
    int fq_poly_is_irreducible_ben_or(const fq_poly_t f, const fq_ctx_t ctx)
    int _fq_poly_is_squarefree(const fq_struct * f, slong len, const fq_ctx_t ctx)
    int fq_poly_is_squarefree(const fq_poly_t f, const fq_ctx_t ctx)
    int fq_poly_factor_equal_deg_prob(fq_poly_t factor, flint_rand_t state, const fq_poly_t pol, slong d, const fq_ctx_t ctx)
    void fq_poly_factor_equal_deg(fq_poly_factor_t factors, const fq_poly_t pol, slong d, const fq_ctx_t ctx)
    void fq_poly_factor_split_single(fq_poly_t linfactor, const fq_poly_t input, const fq_ctx_t ctx)
    void fq_poly_factor_distinct_deg(fq_poly_factor_t res, const fq_poly_t poly, slong * const * degs, const fq_ctx_t ctx)
    void fq_poly_factor_squarefree(fq_poly_factor_t res, const fq_poly_t f, const fq_ctx_t ctx)
    void fq_poly_factor(fq_poly_factor_t res, fq_t lead, const fq_poly_t f, const fq_ctx_t ctx)
    void fq_poly_factor_cantor_zassenhaus(fq_poly_factor_t res, const fq_poly_t f, const fq_ctx_t ctx)
    void fq_poly_factor_kaltofen_shoup(fq_poly_factor_t res, const fq_poly_t poly, const fq_ctx_t ctx)
    void fq_poly_factor_berlekamp(fq_poly_factor_t factors, const fq_poly_t f, const fq_ctx_t ctx)
    void fq_poly_factor_with_berlekamp(fq_poly_factor_t res, fq_t leading_coeff, const fq_poly_t f, const fq_ctx_t ctx)
    void fq_poly_factor_with_cantor_zassenhaus(fq_poly_factor_t res, fq_t leading_coeff, const fq_poly_t f, const fq_ctx_t ctx)
    void fq_poly_factor_with_kaltofen_shoup(fq_poly_factor_t res, fq_t leading_coeff, const fq_poly_t f, const fq_ctx_t ctx)
    void fq_poly_iterated_frobenius_preinv(fq_poly_t * rop, slong n, const fq_poly_t v, const fq_poly_t vinv, const fq_ctx_t ctx)
    void fq_poly_roots(fq_poly_factor_t r, const fq_poly_t f, int with_multiplicity, const fq_ctx_t ctx)
