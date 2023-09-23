from flint.flintlib.arb cimport arb_t
from flint.flintlib.acb cimport acb_t, acb_ptr
from flint.flintlib.fmpz_poly cimport fmpz_poly_t
from flint.flintlib.flint cimport ulong, slong
from flint.flintlib.fmpz cimport fmpz_struct

cdef extern from "flint/arb_fmpz_poly.h":
# from here on is parsed
    void _arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz_struct * poly, slong len, const arb_t x, slong prec)
    void arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)
    void _arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz_struct * poly, slong len, const arb_t x, slong prec)
    void arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)
    void _arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_struct * poly, slong len, const arb_t x, slong prec)
    void arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec)
    void _arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz_struct * poly, slong len, const acb_t x, slong prec)
    void arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)
    void _arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz_struct * poly, slong len, const acb_t x, slong prec)
    void arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)
    void _arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_struct * poly, slong len, const acb_t x, slong prec)
    void arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t poly, const acb_t x, slong prec)
    ulong arb_fmpz_poly_deflation(const fmpz_poly_t poly)
    void arb_fmpz_poly_deflate(fmpz_poly_t res, const fmpz_poly_t poly, ulong deflation)
    void arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, slong prec)
    void arb_fmpz_poly_cos_minpoly(fmpz_poly_t res, ulong n)
    void arb_fmpz_poly_gauss_period_minpoly(fmpz_poly_t res, ulong q, ulong n)
