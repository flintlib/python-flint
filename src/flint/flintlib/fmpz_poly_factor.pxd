from flint.flintlib.fmpz_poly cimport fmpz_poly_t, fmpz_poly_factor_t

cdef extern from "flint/fmpz_poly_factor.h":
    void fmpz_poly_factor_init(fmpz_poly_factor_t fac)
    void fmpz_poly_factor_clear(fmpz_poly_factor_t fac)
    void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t fac, fmpz_poly_t G)
    void fmpz_poly_factor(fmpz_poly_factor_t fac, fmpz_poly_t G)
    void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, fmpz_poly_t G)
