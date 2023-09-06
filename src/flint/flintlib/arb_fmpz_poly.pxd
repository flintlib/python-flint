from flint.flintlib.arb cimport arb_t
from flint.flintlib.acb cimport acb_t, acb_ptr
from flint.flintlib.fmpz_poly cimport fmpz_poly_t
from flint._flint cimport ulong

cdef extern from "arb_fmpz_poly.h":
    void arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t poly, const arb_t x, long prec)
    void arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t poly, const acb_t x, long prec)
    void arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, long prec)
    ulong arb_fmpz_poly_deflation(const fmpz_poly_t poly)
    void arb_fmpz_poly_deflate(fmpz_poly_t res, const fmpz_poly_t poly, ulong deflation)
