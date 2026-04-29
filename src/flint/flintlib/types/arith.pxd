from flint.flintlib.types.flint cimport ulong
from flint.flintlib.types.fmpz cimport fmpz_poly_t
from flint.flintlib.functions.fmpq_poly cimport fmpq_poly_t


cdef extern from "flint/arith.h":
    # Macros
    void arith_legendre_polynomial(fmpq_poly_t v, ulong n)
    void arith_chebyshev_t_polynomial(fmpz_poly_t v, ulong n)
    void arith_chebyshev_u_polynomial(fmpz_poly_t v, ulong n)
    void arith_cyclotomic_polynomial(fmpz_poly_t v, ulong n)
