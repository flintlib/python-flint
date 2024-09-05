from flint.flintlib.types.flint cimport fmpz_t, slong
from flint.flintlib.types.fmpz cimport fmpz_poly_factor_t, fmpz_poly_t



cdef extern from "flint/fmpz_poly_factor.h":
    void fmpz_poly_factor_init(fmpz_poly_factor_t fac)
    void fmpz_poly_factor_init2(fmpz_poly_factor_t fac, slong alloc)
    void fmpz_poly_factor_realloc(fmpz_poly_factor_t fac, slong alloc)
    void fmpz_poly_factor_fit_length(fmpz_poly_factor_t fac, slong len)
    void fmpz_poly_factor_clear(fmpz_poly_factor_t fac)
    void fmpz_poly_factor_set(fmpz_poly_factor_t res, const fmpz_poly_factor_t fac)
    void fmpz_poly_factor_insert(fmpz_poly_factor_t fac, const fmpz_poly_t p, slong e)
    void fmpz_poly_factor_concat(fmpz_poly_factor_t res, const fmpz_poly_factor_t fac)
    void fmpz_poly_factor_print(const fmpz_poly_factor_t fac)
    void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, const fmpz_poly_t F)
    void fmpz_poly_factor_zassenhaus_recombination(fmpz_poly_factor_t final_fac, const fmpz_poly_factor_t lifted_fac, const fmpz_poly_t F, const fmpz_t P, slong exp)
    void _fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, slong exp, const fmpz_poly_t f, slong cutoff, int use_van_hoeij)
    void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, const fmpz_poly_t F)
    void _fmpz_poly_factor_quadratic(fmpz_poly_factor_t fac, const fmpz_poly_t f, slong exp)
    void _fmpz_poly_factor_cubic(fmpz_poly_factor_t fac, const fmpz_poly_t f, slong exp)
    void fmpz_poly_factor(fmpz_poly_factor_t final_fac, const fmpz_poly_t F)
