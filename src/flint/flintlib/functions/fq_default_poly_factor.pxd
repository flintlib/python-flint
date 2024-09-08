from flint.flintlib.types.flint cimport slong, ulong
from flint.flintlib.types.fq_default cimport fq_default_ctx_t, fq_default_poly_factor_t, fq_default_poly_t, fq_default_t



cdef extern from "flint/fq_default_poly_factor.h":
    void fq_default_poly_factor_init(fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_clear(fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_realloc(fq_default_poly_factor_t fac, slong alloc, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_fit_length(fq_default_poly_factor_t fac, slong len, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_set(fq_default_poly_factor_t res, const fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_print_pretty(const fq_default_poly_factor_t fac, const char * var, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_print(const fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_insert(fq_default_poly_factor_t fac, const fq_default_poly_t poly, slong exp, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_concat(fq_default_poly_factor_t res, const fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_pow(fq_default_poly_factor_t fac, slong exp, const fq_default_ctx_t ctx)
    ulong fq_default_poly_remove(fq_default_poly_t f, const fq_default_poly_t p, const fq_default_ctx_t ctx)
    slong fq_default_poly_factor_length(fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_get_poly(fq_default_poly_t poly, const fq_default_poly_factor_t fac, slong i, const fq_default_ctx_t ctx)
    slong fq_default_poly_factor_exp(fq_default_poly_factor_t fac, slong i, const fq_default_ctx_t ctx)
    int fq_default_poly_is_irreducible(const fq_default_poly_t f, const fq_default_ctx_t ctx)
    int fq_default_poly_is_squarefree(const fq_default_poly_t f, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_equal_deg(fq_default_poly_factor_t factors, const fq_default_poly_t pol, slong d, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_split_single(fq_default_poly_t linfactor, const fq_default_poly_t input, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_distinct_deg(fq_default_poly_factor_t res, const fq_default_poly_t poly, slong * const * degs, const fq_default_ctx_t ctx)
    void fq_default_poly_factor_squarefree(fq_default_poly_factor_t res, const fq_default_poly_t f, const fq_default_ctx_t ctx)
    void fq_default_poly_factor(fq_default_poly_factor_t res, fq_default_t lead, const fq_default_poly_t f, const fq_default_ctx_t ctx)
    void fq_default_poly_roots(fq_default_poly_factor_t r, const fq_default_poly_t f, int with_multiplicity, const fq_default_ctx_t ctx)
