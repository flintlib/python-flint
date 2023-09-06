from flint._flint cimport fmpz_struct, ulong
from flint.flintlib.fmpz cimport fmpz_t

cdef extern from "flint/fmpz_factor.h":
    ctypedef struct fmpz_factor_struct:
        int sign
        fmpz_struct * p
        fmpz_struct * exp
        long alloc
        long num
    ctypedef fmpz_factor_struct fmpz_factor_t[1]
    void fmpz_factor_init(fmpz_factor_t factor)
    void fmpz_factor_clear(fmpz_factor_t factor)
    void fmpz_factor(fmpz_factor_t factor, fmpz_t n)
    int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n, ulong start, ulong num_primes)
    void fmpz_factor_expand(fmpz_t n, const fmpz_factor_t factor)
    void _fmpz_factor_append(fmpz_factor_t factor, const fmpz_t p, ulong exp)
    void fmpz_factor_smooth(fmpz_factor_t factor, fmpz_t n, long bits, int proved)
