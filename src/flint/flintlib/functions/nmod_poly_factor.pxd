from flint.flintlib.types.flint cimport flint_rand_t, nmod_t, nn_srcptr, slong, ulong
from flint.flintlib.types.nmod cimport nmod_poly_factor_t, nmod_poly_t



cdef extern from "flint/nmod_poly_factor.h":
    void nmod_poly_factor_init(nmod_poly_factor_t fac)
    void nmod_poly_factor_clear(nmod_poly_factor_t fac)
    void nmod_poly_factor_realloc(nmod_poly_factor_t fac, slong alloc)
    void nmod_poly_factor_fit_length(nmod_poly_factor_t fac, slong len)
    void nmod_poly_factor_set(nmod_poly_factor_t res, const nmod_poly_factor_t fac)
    void nmod_poly_factor_print(const nmod_poly_factor_t fac)
    void nmod_poly_factor_insert(nmod_poly_factor_t fac, const nmod_poly_t poly, slong exp)
    void nmod_poly_factor_concat(nmod_poly_factor_t res, const nmod_poly_factor_t fac)
    void nmod_poly_factor_pow(nmod_poly_factor_t fac, slong exp)
    int nmod_poly_is_irreducible(const nmod_poly_t f)
    int nmod_poly_is_irreducible_ddf(const nmod_poly_t f)
    int nmod_poly_is_irreducible_rabin(const nmod_poly_t f)
    int _nmod_poly_is_squarefree(nn_srcptr f, slong len, nmod_t mod)
    int nmod_poly_is_squarefree(const nmod_poly_t f)
    void nmod_poly_factor_squarefree(nmod_poly_factor_t res, const nmod_poly_t f)
    int nmod_poly_factor_equal_deg_prob(nmod_poly_t factor, flint_rand_t state, const nmod_poly_t pol, slong d)
    void nmod_poly_factor_equal_deg(nmod_poly_factor_t factors, const nmod_poly_t pol, slong d)
    void nmod_poly_factor_distinct_deg(nmod_poly_factor_t res, const nmod_poly_t poly, slong * const * degs)
    void nmod_poly_factor_distinct_deg_threaded(nmod_poly_factor_t res, const nmod_poly_t poly, slong * const * degs)
    void nmod_poly_factor_cantor_zassenhaus(nmod_poly_factor_t res, const nmod_poly_t f)
    void nmod_poly_factor_berlekamp(nmod_poly_factor_t res, const nmod_poly_t f)
    void nmod_poly_factor_kaltofen_shoup(nmod_poly_factor_t res, const nmod_poly_t poly)
    ulong nmod_poly_factor_with_berlekamp(nmod_poly_factor_t res, const nmod_poly_t f)
    ulong nmod_poly_factor_with_cantor_zassenhaus(nmod_poly_factor_t res, const nmod_poly_t f)
    ulong nmod_poly_factor_with_kaltofen_shoup(nmod_poly_factor_t res, const nmod_poly_t f)
    ulong nmod_poly_factor(nmod_poly_factor_t res, const nmod_poly_t f)
    void _nmod_poly_interval_poly_worker(void * arg_ptr)
