from flint._flint cimport ulong
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.fmpz_poly cimport fmpz_poly_t
from flint.flintlib.fmpq_poly cimport fmpq_poly_t
from flint.flintlib.fmpq cimport fmpq_t

cdef extern from "flint/arith.h":
    void arith_number_of_partitions(fmpz_t res, ulong n)
    int arith_moebius_mu(fmpz_t n)
    void arith_divisor_sigma(fmpz_t v, ulong k, fmpz_t n)
    void arith_euler_phi(fmpz_t v, fmpz_t n)
    void arith_bell_number(fmpz_t v, ulong n)
    void arith_euler_number(fmpz_t v, ulong n)
    void arith_bernoulli_number(fmpq_t v, ulong n)
    void arith_stirling_number_1(fmpz_t s, long n, long k)
    void arith_stirling_number_2(fmpz_t s, long n, long k)
    void arith_harmonic_number(fmpq_t v, ulong n)
    void arith_bernoulli_polynomial(fmpq_poly_t v, ulong n)
    void arith_euler_polynomial(fmpq_poly_t v, ulong n)
    void arith_legendre_polynomial(fmpq_poly_t v, ulong n)
    void arith_chebyshev_t_polynomial(fmpz_poly_t v, ulong n)
    void arith_chebyshev_u_polynomial(fmpz_poly_t v, ulong n)
    void arith_cyclotomic_polynomial(fmpz_poly_t v, ulong n)
